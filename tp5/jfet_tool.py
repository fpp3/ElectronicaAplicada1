#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JFET TP Tool: registro I_D–V_DS, detección de IDSS, cálculo de punto Q y rectas de carga.
Stack: tkinter (GUI), csv (IO), math, matplotlib (gráficas estándar).
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog, simpledialog
import tkinter.font as tkfont
import csv, math
from dataclasses import dataclass, field

# Matplotlib embed
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


# -------------------------- Utilidades --------------------------

E_SERIES = {
    "NONE": [],
    "E6":  [10, 15, 22, 33, 47, 68],
    "E12": [10, 12, 15, 18, 22, 27, 33, 39, 47, 56, 68, 82],
    "E24": [10, 11, 12, 13, 15, 16, 18, 20, 22, 24, 27, 30, 33, 36, 39, 43, 47, 51, 56, 62, 68, 75, 82, 91],
}

def normalize_resistance(value_ohm: float, series_name: str) -> float:
    """Redondea value_ohm a la serie E elegida manteniendo década."""
    if series_name == "NONE" or value_ohm <= 0:
        return value_ohm
    series = E_SERIES[series_name]
    decade = 0
    v = value_ohm
    while v >= 100:
        v /= 10
        decade += 1
    while v < 10:
        v *= 10
        decade -= 1
    # elegir el más cercano en la serie
    best = min(series, key=lambda s: abs(s - v))
    return best * (10 ** decade)

def parse_float(entry: tk.Entry, default=None):
    try:
        return float(entry.get().strip().replace(",", "."))
    except:
        return default

def fmt(v, prec=4):
    try:
        return f"{v:.{prec}g}"
    except:
        return str(v)

# -------------------------- Modelo de datos --------------------------

@dataclass
class MeasurementDB:
    samples: list = field(default_factory=list)  # [(VDS, ID)]
    idss: float | None = None                  # en A

    def add_sample(self, vds: float, id_a: float, detect_pct=0.02):
        """Agrega muestra y actualiza IDSS si el cambio relativo < 2%."""
        self.samples.append((vds, id_a))
        if len(self.samples) >= 2:
            prev = self.samples[-2][1]
            if prev > 0:
                rel = abs(id_a - prev) / prev
                if rel < detect_pct:
                    self.idss = id_a
        else:
            # primera muestra no permite detectar IDSS
            pass

    def save_csv(self, path):
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["VDS[V]", "ID[A]"])
            for v, i in self.samples:
                w.writerow([v, i])

    def load_csv(self, path):
        """Carga muestras desde CSV. Soporta encabezados comunes.
        Columnas esperadas (flexible, case-insensitive):
        - VDS: "VDS[V]", "V_DS", "VDS", "VDS (V)", etc.
        - ID:  "ID[A]", "I_D", "ID", "I_D [mA]", etc.
        Si la columna de ID indica mA en el encabezado, convierte a A.
        """
        self.samples.clear()
        self.idss = None

        def _to_float(s):
            try:
                return float(str(s).strip().replace(",", "."))
            except Exception:
                return None

        with open(path, newline="") as f:
            # Intentar detectar dialecto
            try:
                sniff = csv.Sniffer().sniff(f.read(1024))
                f.seek(0)
                reader = csv.reader(f, sniff)
            except Exception:
                f.seek(0)
                reader = csv.reader(f)

            rows = list(reader)
            if not rows:
                return

            header = rows[0]
            data_rows = rows[1:] if header and any(h.strip() for h in header) else rows

            vds_idx = None
            id_idx = None
            id_in_mA = False

            if header and any(h.strip() for h in header):
                hlow = [h.strip().lower() for h in header]
                # Detectar VDS
                for i, h in enumerate(hlow):
                    if ("vds" in h) or ("v_ds" in h) or (h.startswith("v") and "[v]" in h):
                        vds_idx = i
                        break
                # Detectar ID
                for i, h in enumerate(hlow):
                    if (h == "id" or "i_d" in h or "id[" in h or "id " in h):
                        id_idx = i
                        if "ma" in h:
                            id_in_mA = True
                        break

            # Fallback: primeras dos columnas
            if vds_idx is None:
                vds_idx = 0
            if id_idx is None:
                id_idx = 1 if len(data_rows[0]) > 1 else 0

            for row in data_rows:
                if not row or len(row) <= max(vds_idx, id_idx):
                    continue
                v = _to_float(row[vds_idx])
                i = _to_float(row[id_idx])
                if v is None or i is None:
                    continue
                i_a = i / 1000.0 if id_in_mA else i
                self.add_sample(v, i_a)

    def clear(self):
        """Limpia todas las muestras e IDSS."""
        self.samples.clear()
        self.idss = None

    def recompute_idss(self, detect_pct=0.02):
        """Recalcula IDSS usando el criterio de variación relativa."""
        self.idss = None
        if len(self.samples) < 2:
            return
        for idx in range(1, len(self.samples)):
            prev = self.samples[idx-1][1]
            cur = self.samples[idx][1]
            if prev > 0:
                rel = abs(cur - prev) / prev
                if rel < detect_pct:
                    self.idss = cur

# -------------------------- Lógica JFET --------------------------

def shockley_vgs_from_id(id_a, idss_a, vp_v):
    """
    Shockley: ID = IDSS * (1 - VGS/VP)^2
    Devuelve VGS para ID dado. Para n-JFET: VP < 0 => VGS resultará negativo.
    """
    if idss_a is None or idss_a <= 0 or id_a < 0:
        return None
    ratio = id_a / idss_a
    if ratio < 0 or ratio > 1 + 1e-9:
        return None
    root = math.sqrt(max(0.0, ratio))
    vgs = vp_v * (1 - root)  # si VP<0 => VGS<0
    return vgs

def compute_bias_from_idss_vp(idss_a, vp_v, vdd_v, vdsq_v, idq_a=None, series="NONE"):
    """
    Caso autosesgo con RS y RD. Asume V_G=0 (gate a masa por RG grande).
    - Si no dan IDQ, usa IDQ = IDSS/2.
    - Calcula VGSQ por Shockley.
    - RS = -VGSQ/IDQ
    - RD = (VDD - VDSQ)/IDQ - RS
    Retorna dict con valores normalizados y exactos.
    """
    if idss_a is None or idss_a <= 0:
        raise ValueError("IDSS inválido")
    if idq_a is None:
        idq_a = idss_a / 2.0
    vgsq = shockley_vgs_from_id(idq_a, idss_a, vp_v)
    if vgsq is None:
        raise ValueError("IDQ fuera de rango respecto a IDSS/VP")

    rs = -vgsq / idq_a
    rd = (vdd_v - vdsq_v) / idq_a - rs

    gm = None
    if vp_v is not None and vp_v != 0:
        # gm = dI_D/dV_GS evaluado en Q. Shockley => gm = 2*sqrt(IDSS*IDQ)/|VP|
        gm = 2.0 * math.sqrt(max(idss_a, 0.0) * max(idq_a, 0.0)) / abs(vp_v)

    rs_n = normalize_resistance(rs, series)
    rd_n = normalize_resistance(rd, series)
    return {
        "IDQ": idq_a,
        "VGSQ": vgsq,
        "RS": rs, "RD": rd,
        "RS_norm": rs_n, "RD_norm": rd_n,
        "gm": gm,
    }

def loadline_points(vdd_v, rd_ohm, rs_ohm, n=2):
    """
    Recta de carga en el plano (VDS, ID) considerando RS en serie.
    Ecuación DC: VDD = ID*(RD+RS) + VDS  => VDS = VDD - ID*(RD+RS)
    Devuelve 2 puntos extremos para graficar.
    """
    if rd_ohm <= 0 or rs_ohm < 0:
        return [(0, 0), (0, 0)]
    id_intercept = vdd_v / (rd_ohm + rs_ohm)  # VDS=0
    return [(0.0, id_intercept), (vdd_v, 0.0)]

# -------------------------- GUI --------------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("JFET TP Tool")
        self.geometry("1100x700")
        self.db = MeasurementDB()

        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True)

        self.tab_meas = ttk.Frame(self.notebook)
        self.tab_q = ttk.Frame(self.notebook)
        self.tab_load = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_meas, text="1) Mediciones I_D vs V_DS")
        self.notebook.add(self.tab_q, text="2) Punto Q y Resistencias")
        self.notebook.add(self.tab_load, text="3) Recta de Carga")

        self.build_tab_measurements()
        self.build_tab_qpoint()
        self.build_tab_loadline()

    # ---------- Helpers UI ----------
    def make_var_label(self, parent, base: str, sub: str | None = None, unit: str | None = None):
        """Crea un widget tipo etiqueta con base grande y subíndice pequeño.
        Ej.: base='V', sub='DS', unit='[V]' -> V_DS visual con DS reducido.
        """
        frm = ttk.Frame(parent)
        base_font = tkfont.nametofont("TkDefaultFont")
        small_font = tkfont.Font(**base_font.actual())
        small_font.configure(size=max(8, base_font['size'] - 2))

        lbl_base = ttk.Label(frm, text=base, font=base_font)
        lbl_base.pack(side="left", anchor="s")
        if sub:
            lbl_sub = ttk.Label(frm, text=sub, font=small_font)
            lbl_sub.pack(side="left", anchor="s", pady=(6, 0))
        if unit:
            ttk.Label(frm, text=f" {unit}", font=base_font).pack(side="left", anchor="s")
        return frm

    # ---------- Tab 1: Mediciones ----------
    def build_tab_measurements(self):
        fr = self.tab_meas

        # Entradas
        lf_in = ttk.LabelFrame(fr, text="Nueva muestra")
        lf_in.pack(side="left", fill="y", padx=8, pady=8)

        self.make_var_label(lf_in, base="V", sub="DS", unit="[V]").grid(row=0, column=0, sticky="e", padx=4, pady=4)
        self.e_vds = ttk.Entry(lf_in, width=12)
        self.e_vds.grid(row=0, column=1, padx=4, pady=4)

        self.make_var_label(lf_in, base="I", sub="D", unit="[mA]").grid(row=1, column=0, sticky="e", padx=4, pady=4)
        self.e_id_mA = ttk.Entry(lf_in, width=12)
        self.e_id_mA.grid(row=1, column=1, padx=4, pady=4)

        self.btn_add = ttk.Button(lf_in, text="Agregar muestra", command=self.add_sample)
        self.btn_add.grid(row=2, column=0, columnspan=2, pady=6)

        self.btn_load = ttk.Button(lf_in, text="Cargar CSV", command=self.load_csv_ui)
        self.btn_load.grid(row=3, column=0, columnspan=2, pady=6)

        self.btn_save = ttk.Button(lf_in, text="Guardar CSV", command=self.save_csv)
        self.btn_save.grid(row=4, column=0, columnspan=2, pady=6)

        self.btn_clear = ttk.Button(lf_in, text="Limpiar", command=self.clear_samples_ui)
        self.btn_clear.grid(row=5, column=0, columnspan=2, pady=6)

        # Preset de V_DS
        lf_preset = ttk.LabelFrame(lf_in, text="Preset V_DS")
        lf_preset.grid(row=6, column=0, columnspan=2, sticky="ew", padx=2, pady=4)
        self.make_var_label(lf_preset, base="V", sub="DS", unit="inicial [V]").grid(row=0, column=0, sticky="e", padx=4, pady=2)
        self.e_preset_vstart = ttk.Entry(lf_preset, width=10)
        self.e_preset_vstart.insert(0, "0")
        self.e_preset_vstart.grid(row=0, column=1, padx=4, pady=2)

        # ΔV_DS
        frm_delta = ttk.Frame(lf_preset)
        ttk.Label(frm_delta, text="Δ").pack(side="left")
        self.make_var_label(frm_delta, base="V", sub="DS", unit="[V]").pack(side="left")
        frm_delta.grid(row=1, column=0, sticky="e", padx=4, pady=2)
        self.e_preset_step = ttk.Entry(lf_preset, width=10)
        self.e_preset_step.insert(0, "1")
        self.e_preset_step.grid(row=1, column=1, padx=4, pady=2)

        ttk.Label(lf_preset, text="Cantidad").grid(row=2, column=0, sticky="e", padx=4, pady=2)
        self.e_preset_n = ttk.Entry(lf_preset, width=10)
        self.e_preset_n.insert(0, "10")
        self.e_preset_n.grid(row=2, column=1, padx=4, pady=2)

        fr_pbtn = ttk.Frame(lf_preset)
        fr_pbtn.grid(row=3, column=0, columnspan=2, pady=4)
        self.btn_preset_start = ttk.Button(fr_pbtn, text="Iniciar preset", command=self.start_preset)
        self.btn_preset_start.pack(side="left", padx=2)
        self.btn_preset_cancel = ttk.Button(fr_pbtn, text="Cancelar", command=self.cancel_preset)
        self.btn_preset_cancel.pack(side="left", padx=2)

        self.var_preset_status = tk.StringVar(value="")
        ttk.Label(lf_preset, textvariable=self.var_preset_status).grid(row=4, column=0, columnspan=2, sticky="w", padx=4)

        ttk.Separator(lf_in, orient="horizontal").grid(row=7, column=0, columnspan=2, sticky="ew", pady=6)
        ttk.Label(lf_in, text="I_DSS detectado [mA]:").grid(row=8, column=0, sticky="e")
        self.var_idss = tk.StringVar(value="—")
        ttk.Label(lf_in, textvariable=self.var_idss).grid(row=8, column=1, sticky="w")

        ttk.Label(lf_in, text="Tolerancia I_DSS [%]").grid(row=9, column=0, sticky="e", padx=4)
        self.var_detect_pct = tk.StringVar(value="2.0")
        self.spn_detect = tk.Spinbox(lf_in, from_=0.1, to=10.0, increment=0.1, width=6, textvariable=self.var_detect_pct)
        self.spn_detect.grid(row=9, column=1, sticky="w")
        ttk.Button(lf_in, text="Recalcular I_DSS", command=self.recalc_idss_ui).grid(row=10, column=0, columnspan=2, pady=4)

        # Tabla
        lf_tbl = ttk.LabelFrame(fr, text="Muestras registradas")
        lf_tbl.pack(side="top", fill="both", expand=True, padx=8, pady=8)

        self.tree = ttk.Treeview(lf_tbl, columns=("vds", "id"), show="headings", height=10)
        self.tree.heading("vds", text="V_DS [V]")
        self.tree.heading("id", text="I_D [mA]")
        self.tree.column("vds", width=100, anchor="center")
        self.tree.column("id", width=100, anchor="center")
        self.tree.pack(side="left", fill="both", expand=True)
        self.tree.bind("<Double-1>", lambda e: self.edit_selected())

        sb = ttk.Scrollbar(lf_tbl, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=sb.set)
        sb.pack(side="right", fill="y")

        # Controles de edición
        fr_edit = ttk.Frame(lf_tbl)
        fr_edit.pack(fill="x", padx=4, pady=4)
        ttk.Button(fr_edit, text="Editar selección", command=self.edit_selected).pack(side="left", padx=2)
        ttk.Button(fr_edit, text="Eliminar selección", command=self.delete_selected).pack(side="left", padx=2)
        ttk.Button(fr_edit, text="Ordenar por V_DS", command=self.sort_by_vds).pack(side="left", padx=8)

        # Plot
        lf_plot = ttk.LabelFrame(fr, text="Curva I_D vs V_DS")
        lf_plot.pack(side="bottom", fill="both", expand=True, padx=8, pady=8)

        self.fig1 = Figure(figsize=(5, 3), dpi=100)
        self.ax1 = self.fig1.add_subplot(111)
        self.ax1.set_xlabel(r"$V_{DS}$ [V]")
        self.ax1.set_ylabel(r"$I_D$ [mA]")
        self.ax1.grid(True)
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=lf_plot)
        self.canvas1.get_tk_widget().pack(fill="both", expand=True)
        ttk.Button(lf_plot, text="Guardar gráfica", command=self.save_plot1).pack(anchor="e", padx=6, pady=4)

        # Nota de fórmula
        txt = ("Criterio IDSS: al agregar muestras con V_GS≈0, si ΔI_D/I_D(prev) < 2%, "
               "se asume saturación y se fija IDSS = I_D(actual).")
        ttk.Label(fr, text=txt, wraplength=500).pack(side="left", padx=8, pady=8)

    def add_sample(self):
        vds = parse_float(self.e_vds, None)
        id_mA = parse_float(self.e_id_mA, None)
        if vds is None or id_mA is None or id_mA < 0:
            messagebox.showerror("Error", "Ingresá V_DS [V] e I_D [mA] válidos.")
            return
        self.db.add_sample(vds, id_mA / 1000.0, detect_pct=self.get_detect_pct())  # a Amperios
        self.refresh_table()
        self.advance_preset_if_active()

    def save_csv(self):
        if not self.db.samples:
            messagebox.showwarning("Atención", "No hay muestras.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not path:
            return
        self.db.save_csv(path)
        messagebox.showinfo("OK", f"Guardado: {path}")

    def load_csv_ui(self):
        path = filedialog.askopenfilename(filetypes=[("CSV", "*.csv"), ("Todos", "*.*")])
        if not path:
            return
        try:
            self.db.load_csv(path)
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo cargar el CSV:\n{e}")
            return
        self.refresh_table()
        messagebox.showinfo("OK", f"Cargado: {path}")

    def clear_samples_ui(self):
        # Confirmación
        if not messagebox.askyesno("Confirmar", "¿Seguro que querés limpiar todas las muestras? Esta acción no se puede deshacer."):
            return
        # Limpiar datos
        self.db.clear()
        # Refrescar todo
        self.refresh_table()
        # Limpiar entradas de muestra
        self.e_vds.delete(0, "end")
        self.e_id_mA.delete(0, "end")

    # ---------- Preset helpers ----------
    def start_preset(self):
        v0 = parse_float(self.e_preset_vstart, None)
        step = parse_float(self.e_preset_step, None)
        try:
            n = int(float(self.e_preset_n.get().strip()))
        except Exception:
            n = None
        if v0 is None or step is None or step == 0 or n is None or n <= 0:
            messagebox.showerror("Error", "Preset inválido. Revisá V0, ΔV y cantidad.")
            return
        self._preset = {"v0": v0, "step": step, "n": n, "idx": 0, "active": True}
        self.e_vds.delete(0, "end"); self.e_vds.insert(0, fmt(v0))
        self.e_id_mA.delete(0, "end")
        self.var_preset_status.set(f"Preset activo: 1/{n}")

    def cancel_preset(self):
        if hasattr(self, "_preset"):
            self._preset["active"] = False
        self.var_preset_status.set("")

    def advance_preset_if_active(self):
        p = getattr(self, "_preset", None)
        if not p or not p.get("active"):
            return
        p["idx"] += 1
        if p["idx"] >= p["n"]:
            p["active"] = False
            self.var_preset_status.set("Preset completado")
            messagebox.showinfo("Preset", "Se completaron todas las muestras del preset.")
            return
        v_next = p["v0"] + p["idx"] * p["step"]
        self.e_vds.delete(0, "end"); self.e_vds.insert(0, fmt(v_next))
        self.e_id_mA.delete(0, "end")
        self.var_preset_status.set(f"Preset activo: {p['idx']+1}/{p['n']}")

    def refresh_table(self):
        # Tabla con iids alineadas al índice de la muestra
        for item in self.tree.get_children():
            self.tree.delete(item)
        for idx, (v, i) in enumerate(self.db.samples):
            self.tree.insert("", "end", iid=str(idx), values=(fmt(v, 4), fmt(i * 1000.0, 4)))
        # Recalcular IDSS y refrescar plot
        self.db.recompute_idss(self.get_detect_pct())
        self.var_idss.set(fmt(self.db.idss * 1000.0, 4) if self.db.idss is not None else "—")
        self.update_plot_id_vs_vds()

    def recalc_idss_ui(self):
        self.db.recompute_idss(self.get_detect_pct())
        self.var_idss.set(fmt(self.db.idss * 1000.0, 4) if self.db.idss is not None else "—")
        self.update_plot_id_vs_vds()

    def get_detect_pct(self):
        try:
            pct = float(self.var_detect_pct.get()) / 100.0
            return max(0.0, pct)
        except Exception:
            return 0.02

    def delete_selected(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showwarning("Atención", "Seleccioná una fila para eliminar.")
            return
        idx = int(sel[0])
        if not (0 <= idx < len(self.db.samples)):
            return
        if not messagebox.askyesno("Confirmar", "¿Eliminar la muestra seleccionada?"):
            return
        del self.db.samples[idx]
        self.refresh_table()

    def edit_selected(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showwarning("Atención", "Seleccioná una fila para editar.")
            return
        idx = int(sel[0])
        if not (0 <= idx < len(self.db.samples)):
            return
        v_old, i_old_a = self.db.samples[idx]
        # Pedir nuevos valores (como strings); usar los actuales por defecto
        v_str = simpledialog.askstring("Editar V_DS", "Nuevo V_DS [V]:", initialvalue=fmt(v_old, 6), parent=self)
        if v_str is None:
            return
        i_str = simpledialog.askstring("Editar I_D", "Nuevo I_D [mA]:", initialvalue=fmt(i_old_a * 1000.0, 6), parent=self)
        if i_str is None:
            return
        try:
            v_new = float(v_str.strip().replace(",", "."))
            i_new_mA = float(i_str.strip().replace(",", "."))
        except Exception:
            messagebox.showerror("Error", "Valores inválidos.")
            return
        if i_new_mA < 0:
            messagebox.showerror("Error", "I_D no puede ser negativa.")
            return
        self.db.samples[idx] = (v_new, i_new_mA / 1000.0)
        self.refresh_table()

    def sort_by_vds(self):
        if not self.db.samples:
            return
        try:
            self.db.samples.sort(key=lambda t: t[0])
            self.refresh_table()
        except Exception:
            pass

    def update_plot_id_vs_vds(self):
        self.ax1.clear()
        self.ax1.set_xlabel(r"$V_{DS}$ [V]")
        self.ax1.set_ylabel(r"$I_D$ [mA]")
        self.ax1.grid(True)
        if self.db.samples:
            xs = [v for v, _ in self.db.samples]
            ys = [i * 1000.0 for _, i in self.db.samples]
            self.ax1.plot(xs, ys, marker="o")
        self.canvas1.draw_idle()

    # ---------- Tab 2: Punto Q ----------
    def build_tab_qpoint(self):
        fr = self.tab_q

        lf_in = ttk.LabelFrame(fr, text="Parámetros y cálculo")
        lf_in.pack(side="left", fill="y", padx=8, pady=8)

        row = 0
        self.make_var_label(lf_in, base="V", sub="DD", unit="[V]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_vdd = ttk.Entry(lf_in, width=10)
        self.e_vdd.insert(0, "12")
        self.e_vdd.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.make_var_label(lf_in, base="V", sub="DS", unit="(Q) [V]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_vdsq = ttk.Entry(lf_in, width=10)
        self.e_vdsq.insert(0, "6")
        self.e_vdsq.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.make_var_label(lf_in, base="V", sub="P", unit="[V]  (n-JFET es negativo)").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_vp = ttk.Entry(lf_in, width=10)
        self.e_vp.insert(0, "-4")
        self.e_vp.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.make_var_label(lf_in, base="I", sub="DSS", unit="[mA]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_idss_mA = ttk.Entry(lf_in, width=10)
        self.e_idss_mA.insert(0, "")
        self.e_idss_mA.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.make_var_label(lf_in, base="I", sub="D", unit="(Q) [mA]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_idq_mA = ttk.Entry(lf_in, width=10)
        self.e_idq_mA.insert(0, "")
        self.e_idq_mA.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        ttk.Label(lf_in, text="Serie normalizada").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.cmb_series = ttk.Combobox(lf_in, values=list(E_SERIES.keys()), width=8, state="readonly")
        self.cmb_series.set("E24")
        self.cmb_series.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.btn_calc_q = ttk.Button(lf_in, text="Calcular punto Q y resistencias", command=self.calc_qpoint)
        self.btn_calc_q.grid(row=row, column=0, columnspan=2, pady=6)

        # Resultados
        lf_out = ttk.LabelFrame(fr, text="Resultados")
        lf_out.pack(side="left", fill="both", expand=True, padx=8, pady=8)

        self.txt_q = tk.Text(lf_out, height=16)
        self.txt_q.pack(fill="both", expand=True, padx=4, pady=4)

        # Fórmulas usadas
        lf_fx = ttk.LabelFrame(fr, text="Fórmulas usadas")
        lf_fx.pack(side="bottom", fill="x", padx=8, pady=8)
        formulas = (
            "Shockley:  I_D = I_DSS · (1 - V_GS/VP)^2\n"
            "⇒ V_GS(Q) = VP · (1 - √(I_D(Q)/I_DSS))\n"
            "Autosesgo (V_G≈0):  V_GS(Q) = -I_D(Q)·R_S\n"
            "Lazo DC:  V_DD = I_D(Q)·(R_D + R_S) + V_DS(Q)\n"
            "⇒ R_S = -V_GS(Q)/I_D(Q)\n"
            "⇒ R_D = (V_DD - V_DS(Q))/I_D(Q) - R_S"
        )
        ttk.Label(lf_fx, text=formulas, justify="left").pack(anchor="w", padx=6, pady=6)

    def calc_qpoint(self):
        vdd = parse_float(self.e_vdd, None)
        vdsq = parse_float(self.e_vdsq, None)
        vp = parse_float(self.e_vp, None)
        idss_mA = parse_float(self.e_idss_mA, None)

        # Si no se ingresó IDSS, usar el detectado en Tab 1
        if idss_mA is None and self.db.idss is not None:
            idss_a = self.db.idss
        elif idss_mA is not None:
            idss_a = idss_mA / 1000.0
        else:
            messagebox.showerror("Error", "Ingresá IDSS [mA] o detectalo en la pestaña 1.")
            return

        idq_mA = parse_float(self.e_idq_mA, None)
        idq_a = idq_mA / 1000.0 if idq_mA is not None and idq_mA > 0 else None

        if None in (vdd, vdsq, vp):
            messagebox.showerror("Error", "Completá VDD, VDSQ y VP.")
            return

        series = self.cmb_series.get()
        try:
            r = compute_bias_from_idss_vp(idss_a, vp, vdd, vdsq, idq_a, series)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        rg = parse_float(self.e_rg, None) if hasattr(self, "e_rg") else None
        if rg is None or rg <= 0:
            zin_val = math.inf
            zin_label = "~inf (sin RG)"
        else:
            zin_val = rg
            zin_label = f"{fmt(zin_val)} Ω (domina RG)"

        rd_used = r['RD_norm'] if series != "NONE" and r['RD_norm'] > 0 else r['RD']
        zout_val = rd_used
        zout_label = f"{fmt(zout_val)} Ω (domina RD)"

        # Salida
        self.txt_q.delete("1.0", "end")
        lines = []
        lines.append(f"I_DSS   = {fmt(idss_a*1000)} mA")
        lines.append(f"I_D(Q)  = {fmt(r['IDQ']*1000)} mA")
        lines.append(f"V_GS(Q) = {fmt(r['VGSQ'])} V")
        lines.append(f"RS   = {fmt(r['RS'])} Ω   → RS_norm({series}) = {fmt(r['RS_norm'])} Ω")
        lines.append(f"RD   = {fmt(r['RD'])} Ω   → RD_norm({series}) = {fmt(r['RD_norm'])} Ω")
        if r.get("gm") is not None:
            lines.append(f"g_m(Q) = {fmt(r['gm'])} S")
        lines.append(f"Z_entrada approx {zin_label}")
        lines.append(f"Z_salida approx {zout_label}")
        self.txt_q.insert("end", "\n".join(lines))

        # Guardar últimos para Tab 3
        self.last_rd = r['RD_norm'] if series != "NONE" else r['RD']
        self.last_rs = r['RS_norm'] if series != "NONE" else r['RS']
        self.last_idq = r['IDQ']
        self.last_vdsq = vdsq
        self.last_zin = zin_val
        self.last_zout = zout_val
        if hasattr(self, "var_zin"):
            self.var_zin.set(zin_label)
        if hasattr(self, "var_zout"):
            self.var_zout.set(zout_label)
        self.notebook.select(self.tab_load)

    # ---------- Tab 3: Recta de carga ----------
    def build_tab_loadline(self):
        fr = self.tab_load

        lf_in = ttk.LabelFrame(fr, text="Parámetros de la etapa y mediciones")
        lf_in.pack(side="left", fill="y", padx=8, pady=8)

        row = 0
        self.make_var_label(lf_in, base="V", sub="DD", unit="[V]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_vdd2 = ttk.Entry(lf_in, width=10)
        self.e_vdd2.insert(0, "12")
        self.e_vdd2.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.make_var_label(lf_in, base="R", sub="D", unit="[Ω]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_rd = ttk.Entry(lf_in, width=10)
        self.e_rd.insert(0, "")
        self.e_rd.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.make_var_label(lf_in, base="R", sub="S", unit="[Ω]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_rs = ttk.Entry(lf_in, width=10)
        self.e_rs.insert(0, "")
        self.e_rs.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.make_var_label(lf_in, base="R", sub="G", unit="[Ω] (opcional)").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_rg = ttk.Entry(lf_in, width=10)
        self.e_rg.insert(0, "1000000")
        self.e_rg.grid(row=row, column=1, padx=4, pady=4)

        ttk.Separator(lf_in, orient="horizontal").grid(row=row+1, column=0, columnspan=2, sticky="ew", pady=6)

        row += 2
        self.make_var_label(lf_in, base="V", sub="DS", unit="(Q) medido [V]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_vdsq_meas = ttk.Entry(lf_in, width=10); self.e_vdsq_meas.insert(0, "")
        self.e_vdsq_meas.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.make_var_label(lf_in, base="I", sub="D", unit="(Q) medido [mA]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_idq_meas = ttk.Entry(lf_in, width=10); self.e_idq_meas.insert(0, "")
        self.e_idq_meas.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.make_var_label(lf_in, base="V", sub="GS", unit="(Q) medido [V]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.e_vgsq_meas = ttk.Entry(lf_in, width=10); self.e_vgsq_meas.insert(0, "")
        self.e_vgsq_meas.grid(row=row, column=1, padx=4, pady=4)

        row += 1
        self.btn_plot = ttk.Button(lf_in, text="Graficar recta de carga", command=self.plot_loadline)
        self.btn_plot.grid(row=row, column=0, columnspan=2, pady=6)

        row += 1
        ttk.Separator(lf_in, orient="horizontal").grid(row=row, column=0, columnspan=2, sticky="ew", pady=6)

        row += 1
        ttk.Label(lf_in, text="Z_entrada estimada [Ω]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.var_zin = tk.StringVar(value="—")
        ttk.Label(lf_in, textvariable=self.var_zin).grid(row=row, column=1, sticky="w", padx=4, pady=4)

        row += 1
        ttk.Label(lf_in, text="Z_salida estimada [Ω]").grid(row=row, column=0, sticky="e", padx=4, pady=4)
        self.var_zout = tk.StringVar(value="—")
        ttk.Label(lf_in, textvariable=self.var_zout).grid(row=row, column=1, sticky="w", padx=4, pady=4)

        # Plot
        lf_plot = ttk.LabelFrame(fr, text="Recta de carga en (V_DS, I_D)")
        lf_plot.pack(side="left", fill="both", expand=True, padx=8, pady=8)

        self.fig2 = Figure(figsize=(5, 4), dpi=100)
        self.ax2 = self.fig2.add_subplot(111)
        self.ax2.set_xlabel(r"$V_{DS}$ [V]")
        self.ax2.set_ylabel(r"$I_D$ [mA]")
        self.ax2.grid(True)
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=lf_plot)
        self.canvas2.get_tk_widget().pack(fill="both", expand=True)
        ttk.Button(lf_plot, text="Guardar gráfica", command=self.save_plot2).pack(anchor="e", padx=6, pady=4)

        # Tips
        tips = (
            "Recta: V_DS = V_DD - I_D·(R_D+R_S)\n"
            "Puntos usados: (0, V_DD/(R_D+R_S)) y (V_DD, 0)\n"
            "Podés superponer el Punto Q medido u obtenido en la pestaña 2."
        )
        ttk.Label(fr, text=tips, justify="left").pack(side="bottom", anchor="w", padx=8, pady=6)

    def plot_loadline(self):
        vdd = parse_float(self.e_vdd2, None)
        rd = parse_float(self.e_rd, None)
        rs = parse_float(self.e_rs, None)

        # Autocompletar con últimos calculados si están vacíos
        if (rd is None or rd <= 0) and hasattr(self, "last_rd"):
            rd = self.last_rd
            self.e_rd.delete(0, "end"); self.e_rd.insert(0, fmt(rd))
        if (rs is None or rs < 0) and hasattr(self, "last_rs"):
            rs = self.last_rs
            self.e_rs.delete(0, "end"); self.e_rs.insert(0, fmt(rs))

        if None in (vdd, rd, rs) or rd <= 0 or rs < 0:
            messagebox.showerror("Error", "Completá VDD, RD y RS válidos.")
            return

        rg = parse_float(self.e_rg, None)
        if rg is None or rg <= 0:
            zin_val = math.inf
            zin_txt = "~inf (sin RG)"
        else:
            zin_val = rg
            zin_txt = f"{fmt(zin_val)} Ω"

        pts = loadline_points(vdd, rd, rs)

        self.ax2.clear()
        self.ax2.set_xlabel(r"$V_{DS}$ [V]")
        self.ax2.set_ylabel(r"$I_D$ [mA]")
        self.ax2.grid(True)

        xs = [pts[0][0], pts[1][0]]
        ys_mA = [pts[0][1]*1000.0, pts[1][1]*1000.0]
        self.ax2.plot(xs, ys_mA, marker="o", label="Recta de carga")

        # Punto Q calculado previo
        if hasattr(self, "last_idq") and hasattr(self, "last_vdsq"):
            self.ax2.plot(self.last_vdsq, self.last_idq*1000.0, "s", label="Q calc")

        # Punto Q medido
        vdsq_m = parse_float(self.e_vdsq_meas, None)
        idq_m_mA = parse_float(self.e_idq_meas, None)
        if vdsq_m is not None and idq_m_mA is not None:
            self.ax2.plot(vdsq_m, idq_m_mA, "x", label="Q medido")

        self.ax2.legend()
        self.canvas2.draw_idle()

        zout_val = rd
        zout_txt = f"{fmt(zout_val)} Ω"
        self.var_zin.set(zin_txt)
        self.var_zout.set(zout_txt)
        self.last_zin = zin_val
        self.last_zout = zout_val

    def save_plot1(self):
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if not path:
            return
        try:
            self.fig1.savefig(path, dpi=150, bbox_inches="tight")
            messagebox.showinfo("OK", f"Gráfica guardada: {path}")
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo guardar la gráfica:\n{e}")

    def save_plot2(self):
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if not path:
            return
        try:
            self.fig2.savefig(path, dpi=150, bbox_inches="tight")
            messagebox.showinfo("OK", f"Gráfica guardada: {path}")
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo guardar la gráfica:\n{e}")


if __name__ == "__main__":
    App().mainloop()
