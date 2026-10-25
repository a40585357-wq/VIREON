#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
"""
VIREON v3 - Cellular & Synthetic Biology Suite (single file)
Requires: matplotlib, numpy, scipy, rdkit
Install: pip3.8 install matplotlib numpy scipy rdkit
"""

import os
import sys
import math
import threading
import urllib.request
import urllib.parse
import ssl
import json
import random

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from scipy.integrate import odeint

# ─── Optional scientific libraries (Category C) ───
# SPEED: heavy optional libraries (cobra, rdkit, biopython, pandas, ...) are
# NO LONGER imported here at startup. They are lazy-loaded via the _ensure_*()
# helpers defined in biosim_part2 and pulled in through `from biosim_part2
# import *` below. This makes the app start almost instantly instead of taking
# 25-45 seconds loading libraries the language screen never uses.
# The HAS_* flags and _ensure_*() functions all come from biosim_part2.

# Optional Pillow for loading PNG human image (lightweight, keep eager)
try:
    from PIL import Image, ImageTk
    HAS_PIL = True
except ImportError as _exc:
    # biosim_log lives in biosim_part2, which is imported below — not yet
    # available here, so guard the call to avoid a NameError when PIL is absent.
    if 'biosim_log' in globals():
        biosim_log('warn', 'module: {}'.format(_exc))
    HAS_PIL = False


# ─────────────────────────────────────────────────
# LANGUAGES
# ─────────────────────────────────────────────────

# ─── Import shared code + SynBio/MolBio windows from part 2 ───
from biosim_part2 import *
import biosim_part2
from biosim_part3 import *
import biosim_part3

class LangScreen(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('VIREON')
        self.geometry('480x340')
        self.configure(bg=Theme.c('bg'))
        self.resizable(False, False)
        self.lv = tk.StringVar(value='ru')
        self._build()
        self.eval('tk::PlaceWindow . center')
        # Quit cleanly if the start window is closed
        self.protocol('WM_DELETE_WINDOW', self._quit_all)

    def _quit_all(self):
        try:
            self.destroy()
        except Exception as _exc:
            biosim_log('warn', '_quit_all: {}'.format(_exc))
        import os
        os._exit(0)

    def _build(self):
        tk.Label(self, text='VIREON', font=('Arial', 32, 'bold'),
                 bg=Theme.c('bg'), fg='#0969da').pack(pady=(28, 4))
        tk.Label(self, text='Cellular & Synthetic Biology Suite', font=('Arial', 13),
                 bg=Theme.c('bg'), fg=Theme.c('text_dim')).pack()
        # AI status indicators (torch/deepchem load lazily when needed)
        status_frame = tk.Frame(self, bg=Theme.c('bg'))
        status_frame.pack(pady=(4, 8))
        for name, ok in [('RDKit', HAS_RDKIT), ('PIL', HAS_PIL)]:
            txt = '{}: {}'.format(name, 'ON' if ok else 'OFF')
            tk.Label(status_frame, text=txt, font=('Arial', 8),
                     bg=Theme.c('bg'), fg='#1a7f37' if ok else '#57606a').pack(side='left', padx=4)
        tk.Label(self, text='Choose language / Выберите язык', font=('Arial', 12),
                 bg=Theme.c('bg'), fg=Theme.c('text')).pack(pady=(8, 12))
        f = tk.Frame(self, bg=Theme.c('bg')); f.pack()
        for lbl, val in [('Русский', 'ru'), ('English', 'en')]:
            tk.Radiobutton(f, text=lbl, variable=self.lv, value=val,
                           font=('Arial', 13), bg=Theme.c('bg'), fg=Theme.c('text'),
                           selectcolor='#ffffff', activebackground='#f5f7fa',
                           activeforeground='#0969da').pack(anchor='w', pady=4)
        tk.Button(self, text='Start / Запустить', font=('Arial', 14, 'bold'),
                  bg='#1a7f37', highlightbackground='#1a7f37', fg='white', activebackground='#2da44e',
                  relief='flat', cursor='hand2', padx=24, pady=8,
                  command=self._start).pack(pady=(18, 0))

    def _start(self):
        lang = self.lv.get()
        self.destroy()
        MainApp(lang).mainloop()


# ─────────────────────────────────────────────────
# NEW DISEASE DIALOG (extended)
# ─────────────────────────────────────────────────
class NewDiseaseDialog(tk.Toplevel):
    def __init__(self, parent, lang, on_submit):
        super().__init__(parent)
        self.lang = lang
        self.L = LANG[lang]
        self.on_submit = on_submit
        self.title(self.L['add_new_dis'])
        self.geometry('520x600')
        self.configure(bg=Theme.c('bg'))
        self.resizable(False, False)
        self.type_var = tk.StringVar(value='infection')
        self._build()
        self.transient(parent)
        self.grab_set()

    def _build(self):
        tk.Label(self, text=self.L['add_new_dis'],
                 font=('Arial', 14, 'bold'),
                 bg=Theme.c('bg'), fg='#0969da').pack(pady=(12, 8))

        tk.Label(self, text=self.L['dis_name'] + ':',
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                 font=('Arial', 10)).pack(anchor='w', padx=16)
        self.name_var = tk.StringVar()
        tk.Entry(self, textvariable=self.name_var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                 insertbackground='#1c2128', relief='flat', font=('Arial', 11),
                 highlightthickness=1, highlightbackground=Theme.c('border')
                 ).pack(fill='x', padx=16, pady=(2, 6), ipady=4)

        tk.Label(self, text=self.L['dis_type'] + ':',
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                 font=('Arial', 10)).pack(anchor='w', padx=16)
        type_frame = tk.Frame(self, bg=Theme.c('bg'))
        type_frame.pack(fill='x', padx=16, pady=(2, 6))
        for txt, val in [('Инфекция/Infection','infection'),
                          ('Онкология/Cancer','cancer'),
                          ('Генетика/Genetic','genetic'),
                          ('Другое/Other','other')]:
            tk.Radiobutton(type_frame, text=txt, variable=self.type_var, value=val,
                            bg=Theme.c('bg'), fg=Theme.c('text'), selectcolor='#ffffff',
                            activebackground='#f5f7fa', font=('Arial', 9)
                            ).pack(side='left', padx=2)

        tk.Label(self, text=self.L['dis_desc'] + ':',
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                 font=('Arial', 10)).pack(anchor='w', padx=16)
        self.desc_text = tk.Text(self, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                  insertbackground='#1c2128', relief='flat',
                                  font=('Arial', 10), height=3,
                                  highlightthickness=1, highlightbackground=Theme.c('border'))
        self.desc_text.pack(fill='x', padx=16, pady=(2, 6))

        tk.Label(self, text=self.L['dis_genes'] + ':',
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                 font=('Arial', 10)).pack(anchor='w', padx=16)
        tk.Label(self, text='Example: TP53, BRCA1, EGFR',
                 bg=Theme.c('bg'), fg='#484f58',
                 font=('Arial', 8, 'italic')).pack(anchor='w', padx=16)
        self.genes_var = tk.StringVar()
        tk.Entry(self, textvariable=self.genes_var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                 insertbackground='#1c2128', relief='flat', font=('Arial', 11),
                 highlightthickness=1, highlightbackground=Theme.c('border')
                 ).pack(fill='x', padx=16, pady=(2, 6), ipady=4)

        tk.Label(self, text=self.L['dis_symptoms'] + ':',
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                 font=('Arial', 10)).pack(anchor='w', padx=16)
        self.symptoms_var = tk.StringVar()
        tk.Entry(self, textvariable=self.symptoms_var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                 insertbackground='#1c2128', relief='flat', font=('Arial', 11),
                 highlightthickness=1, highlightbackground=Theme.c('border')
                 ).pack(fill='x', padx=16, pady=(2, 6), ipady=4)

        tk.Label(self, text=self.L['dis_organs'] + ':',
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                 font=('Arial', 10)).pack(anchor='w', padx=16)
        tk.Label(self, text='Example: Lungs, Heart, Liver',
                 bg=Theme.c('bg'), fg='#484f58',
                 font=('Arial', 8, 'italic')).pack(anchor='w', padx=16)
        self.organs_var = tk.StringVar()
        tk.Entry(self, textvariable=self.organs_var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                 insertbackground='#1c2128', relief='flat', font=('Arial', 11),
                 highlightthickness=1, highlightbackground=Theme.c('border')
                 ).pack(fill='x', padx=16, pady=(2, 10), ipady=4)

        bf = tk.Frame(self, bg=Theme.c('bg')); bf.pack(fill='x', padx=16, pady=10)
        tk.Button(bf, text=self.L['submit'], bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                  relief='flat', font=('Arial', 11, 'bold'), cursor='hand2',
                  pady=8, command=self._do_submit).pack(side='left', fill='x', expand=True, padx=(0, 4))
        tk.Button(bf, text=self.L['cancel'], bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                  relief='solid', bd=1, font=('Arial', 11), cursor='hand2',
                  pady=8, command=self.destroy).pack(side='right', padx=(4, 0))

    def _do_submit(self):
        name = self.name_var.get().strip()
        desc = self.desc_text.get('1.0', 'end').strip()
        genes = self.genes_var.get().strip()
        symptoms = self.symptoms_var.get().strip()
        organs = self.organs_var.get().strip()
        if not name:
            messagebox.showinfo('VIREON', 'Введите название болезни')
            return
        if not genes:
            messagebox.showinfo('VIREON', 'Введите хотя бы один ген')
            return
        self.destroy()
        self.on_submit(name, desc, genes, symptoms, organs, self.type_var.get())


# ─────────────────────────────────────────────────
# MAIN APP
# ─────────────────────────────────────────────────

def pk_one_compartment(dose_mg, ka, ke, V, time_points, route='oral'):
    """One-compartment PK model.
    
    Args:
        dose_mg: dose in mg
        ka: absorption rate constant (1/h) — for oral
        ke: elimination rate constant (1/h)
        V: volume of distribution (L)
        time_points: array of time values (hours)
        route: 'oral' or 'iv'
    
    Returns:
        Array of concentrations (mg/L) at each time point.
    """
    import math
    concentrations = []
    for t in time_points:
        if route == 'iv':
            # IV bolus: C = (Dose/V) * exp(-ke*t)
            C = (dose_mg / V) * math.exp(-ke * t)
        else:
            # Oral: C = (F*Dose*ka)/(V*(ka-ke)) * (exp(-ke*t) - exp(-ka*t))
            F = 1.0  # bioavailability
            if abs(ka - ke) < 1e-9:
                # Avoid divide by zero
                C = (F * dose_mg / V) * ka * t * math.exp(-ka * t)
            else:
                C = (F * dose_mg * ka) / (V * (ka - ke)) * (
                    math.exp(-ke * t) - math.exp(-ka * t))
        concentrations.append(max(0, C))
    return concentrations

def pk_parameters(time_points, concentrations):
    """Calculate PK parameters from concentration-time data.
    
    Returns dict with:
        Cmax: peak concentration
        Tmax: time to peak
        AUC: area under curve (trapezoid rule)
        t_half: half-life
        Cl: clearance estimate
    """
    if not concentrations or not time_points:
        return {}
    
    # Cmax / Tmax
    cmax_idx = concentrations.index(max(concentrations))
    Cmax = concentrations[cmax_idx]
    Tmax = time_points[cmax_idx]
    
    # AUC (trapezoid)
    AUC = 0
    for i in range(len(time_points) - 1):
        dt = time_points[i+1] - time_points[i]
        AUC += dt * (concentrations[i] + concentrations[i+1]) / 2
    
    # Half-life from terminal phase (last 4 points)
    import math
    t_half = 0
    if len(time_points) >= 4 and concentrations[-1] > 0:
        # ln(C) vs t — slope = -ke
        n = min(4, len(time_points))
        x_vals = time_points[-n:]
        y_vals = [math.log(max(1e-9, c)) for c in concentrations[-n:]]
        # Linear regression slope
        n_pts = len(x_vals)
        sx = sum(x_vals)
        sy = sum(y_vals)
        sxy = sum(x * y for x, y in zip(x_vals, y_vals))
        sxx = sum(x * x for x in x_vals)
        denom = n_pts * sxx - sx * sx
        if abs(denom) > 1e-9:
            slope = (n_pts * sxy - sx * sy) / denom
            if slope < 0:
                t_half = math.log(2) / abs(slope)
    
    return {
        'Cmax': Cmax,
        'Tmax': Tmax,
        'AUC': AUC,
        't_half': t_half,
    }

def pd_emax_model(C, Emax, EC50, n=1):
    """Sigmoid Emax PD model (Hill equation).
    
    Effect = Emax * C^n / (EC50^n + C^n)
    """
    if C <= 0:
        return 0
    return Emax * (C ** n) / ((EC50 ** n) + (C ** n))


# Kaplan-Meier survival

def kaplan_meier(times, events):
    """Compute Kaplan-Meier survival curve.
    
    Args:
        times: list of observation times
        events: list of 1 (event=death) or 0 (censored)
    
    Returns:
        dict with t_points and survival probabilities S(t)
    """
    if not times:
        return {'time': [], 'survival': [], 'at_risk': []}
    
    # Sort by time
    data = sorted(zip(times, events))
    sorted_times = [d[0] for d in data]
    sorted_events = [d[1] for d in data]
    
    unique_times = sorted(set(sorted_times))
    S = 1.0
    t_out = [0]
    s_out = [1.0]
    at_risk_out = [len(sorted_times)]
    
    for t in unique_times:
        # Count at risk at time t
        n_at_risk = sum(1 for tt in sorted_times if tt >= t)
        # Count deaths at time t
        d = sum(1 for tt, e in zip(sorted_times, sorted_events) if tt == t and e == 1)
        if n_at_risk > 0:
            S *= (1 - d / n_at_risk)
        t_out.append(t)
        s_out.append(S)
        at_risk_out.append(n_at_risk - d)
    
    return {
        'time': t_out,
        'survival': s_out,
        'at_risk': at_risk_out,
    }

def log_rank_test(times_a, events_a, times_b, events_b):
    """Log-rank test for comparing two survival curves.
    
    Returns: dict with chi_squared, p_value (approximate), and N per group.
    """
    import math
    # Combine all unique event times
    all_times = sorted(set(times_a + times_b))
    
    O1 = 0  # observed events group 1
    E1 = 0  # expected events group 1
    V = 0   # variance
    
    for t in all_times:
        # At risk
        n1 = sum(1 for tt in times_a if tt >= t)
        n2 = sum(1 for tt in times_b if tt >= t)
        n = n1 + n2
        if n == 0:
            continue
        # Deaths
        d1 = sum(1 for tt, e in zip(times_a, events_a) if tt == t and e == 1)
        d2 = sum(1 for tt, e in zip(times_b, events_b) if tt == t and e == 1)
        d = d1 + d2
        if d == 0:
            continue
        # Expected and variance
        e1 = d * n1 / n
        if n > 1:
            v = (n1 * n2 * d * (n - d)) / (n * n * (n - 1))
        else:
            v = 0
        O1 += d1
        E1 += e1
        V += v
    
    if V < 1e-9:
        return {'chi_squared': 0, 'p_value': 1.0,
                 'n_a': len(times_a), 'n_b': len(times_b)}
    chi_sq = (O1 - E1) ** 2 / V
    # Approximate p-value from chi-square (df=1) using survival function
    # P(X >= x) for chi-square df=1 = 2 * (1 - Φ(sqrt(x)))
    z = math.sqrt(chi_sq)
    # Normal CDF approximation (Abramowitz)
    p = 2 * (1 - 0.5 * (1 + math.erf(z / math.sqrt(2))))
    return {
        'chi_squared': chi_sq,
        'p_value': p,
        'O1': O1, 'E1': E1, 'V': V,
        'n_a': len(times_a), 'n_b': len(times_b),
    }


# Classical statistical tests

def t_test(group_a, group_b, paired=False):
    """Student's t-test (independent or paired).
    
    Returns dict: t, df, p_value, mean_a, mean_b, ci_diff (95%).
    """
    import math
    n_a = len(group_a)
    n_b = len(group_b)
    if n_a < 2 or n_b < 2:
        return None
    mean_a = sum(group_a) / n_a
    mean_b = sum(group_b) / n_b
    var_a = sum((x - mean_a) ** 2 for x in group_a) / (n_a - 1)
    var_b = sum((x - mean_b) ** 2 for x in group_b) / (n_b - 1)
    
    if paired:
        if n_a != n_b:
            return None
        diffs = [a - b for a, b in zip(group_a, group_b)]
        mean_d = sum(diffs) / n_a
        var_d = sum((d - mean_d) ** 2 for d in diffs) / (n_a - 1)
        se = math.sqrt(var_d / n_a)
        t = mean_d / se if se > 0 else 0
        df = n_a - 1
    else:
        # Welch's t-test (unequal variances)
        se = math.sqrt(var_a / n_a + var_b / n_b)
        t = (mean_a - mean_b) / se if se > 0 else 0
        # Welch-Satterthwaite df
        num = (var_a / n_a + var_b / n_b) ** 2
        den = (var_a / n_a) ** 2 / (n_a - 1) + (var_b / n_b) ** 2 / (n_b - 1)
        df = num / den if den > 0 else n_a + n_b - 2
    
    # p-value (two-sided) from t distribution (approximation)
    # Use Student-t survival function
    def t_sf(t, df):
        # Approximation via Beta function: 1 - I_{x}(df/2, 1/2) / 2
        # Using normal approximation for large df
        if df > 30:
            # Normal approximation
            return 1 - 0.5 * (1 + math.erf(abs(t) / math.sqrt(2)))
        # For small df use simpler approximation
        x = df / (df + t * t)
        # Wilson-Hilferty
        return 0.5 * (1 - math.tanh((1 - 2 / (9 * df)) * (1 - x ** (1/3))
                                       / math.sqrt(2 * x ** (1/3) / (9 * df))))
    
    p = 2 * t_sf(t, df)
    p = max(0, min(1, p))
    
    return {
        't': t,
        'df': df,
        'p_value': p,
        'mean_a': mean_a,
        'mean_b': mean_b,
        'std_a': math.sqrt(var_a),
        'std_b': math.sqrt(var_b),
        'n_a': n_a, 'n_b': n_b,
        'diff': mean_a - mean_b,
        'se': se,
    }

def mann_whitney_u(group_a, group_b):
    """Mann-Whitney U test (non-parametric)."""
    import math
    n_a = len(group_a)
    n_b = len(group_b)
    if n_a == 0 or n_b == 0:
        return None
    # Rank all values together
    combined = [(v, 'a') for v in group_a] + [(v, 'b') for v in group_b]
    combined.sort()
    # Assign ranks (handle ties with average)
    ranks = {}
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) - 1 and combined[j+1][0] == combined[i][0]:
            j += 1
        avg_rank = (i + j + 2) / 2  # 1-indexed rank average
        for k in range(i, j + 1):
            ranks[k] = avg_rank
        i = j + 1
    # Sum of ranks for group A
    R_a = sum(ranks[i] for i, (_, grp) in enumerate(combined) if grp == 'a')
    U_a = R_a - n_a * (n_a + 1) / 2
    U_b = n_a * n_b - U_a
    U = min(U_a, U_b)
    # Normal approximation
    mean_U = n_a * n_b / 2
    std_U = math.sqrt(n_a * n_b * (n_a + n_b + 1) / 12)
    if std_U < 1e-9:
        return {'U': U, 'z': 0, 'p_value': 1.0}
    z = (U - mean_U) / std_U
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return {
        'U': U, 'z': z, 'p_value': p,
        'n_a': n_a, 'n_b': n_b,
    }

def anova_oneway(*groups):
    """One-way ANOVA (F-test)."""
    import math
    k = len(groups)
    if k < 2:
        return None
    all_vals = [v for g in groups for v in g]
    N = len(all_vals)
    if N < k + 1:
        return None
    grand_mean = sum(all_vals) / N
    # Between-group SS
    SSB = sum(len(g) * (sum(g) / len(g) - grand_mean) ** 2
               for g in groups if len(g) > 0)
    # Within-group SS
    SSW = 0
    for g in groups:
        if len(g) > 0:
            m = sum(g) / len(g)
            SSW += sum((v - m) ** 2 for v in g)
    df_b = k - 1
    df_w = N - k
    if df_w <= 0 or SSW < 1e-9:
        return None
    MSB = SSB / df_b
    MSW = SSW / df_w
    F = MSB / MSW
    # Approximate p-value
    # F distribution survival function approximation
    # For large df_w, F follows chi-square / df_b
    x = df_w / (df_w + df_b * F)
    p = x ** (df_w / 2)  # very rough approximation
    p = max(0, min(1, p))
    return {
        'F': F, 'df_between': df_b, 'df_within': df_w,
        'p_value': p, 'SSB': SSB, 'SSW': SSW, 'k': k, 'N': N,
    }

def confidence_interval_95(values):
    """95% CI for the mean."""
    import math
    n = len(values)
    if n < 2:
        return None
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    se = math.sqrt(var / n)
    # t_{0.025, df} for various df
    t_crit = 1.96 if n > 30 else 2.0  # rough
    if n == 2: t_crit = 12.71
    elif n == 3: t_crit = 4.30
    elif n == 4: t_crit = 3.18
    elif n == 5: t_crit = 2.78
    elif n <= 10: t_crit = 2.26
    elif n <= 20: t_crit = 2.09
    return {
        'mean': mean,
        'lower': mean - t_crit * se,
        'upper': mean + t_crit * se,
        'se': se,
        'n': n,
    }

def search_pubmed(query, max_results=20, timeout=15):
    """Search PubMed via NCBI E-utilities.
    
    Returns list of dicts: {pmid, title, authors, journal, year, abstract}
    """
    import urllib.request, urllib.parse, json
    try:
        # Step 1: esearch to get PMIDs
        url1 = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
                 '?db=pubmed&term=' + urllib.parse.quote(query) +
                 '&retmax=' + str(max_results) +
                 '&retmode=json')
        with urllib.request.urlopen(url1, timeout=timeout) as f:
            data = json.loads(f.read().decode('utf-8'))
        pmids = data.get('esearchresult', {}).get('idlist', [])
        if not pmids:
            return []
        
        # Step 2: esummary to get titles/authors
        url2 = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
                 '?db=pubmed&id=' + ','.join(pmids) +
                 '&retmode=json')
        with urllib.request.urlopen(url2, timeout=timeout) as f:
            summary = json.loads(f.read().decode('utf-8'))
        
        results = []
        for pmid in pmids:
            entry = summary.get('result', {}).get(pmid, {})
            if not entry:
                continue
            authors = entry.get('authors', [])
            author_names = ', '.join(a.get('name', '') for a in authors[:3])
            if len(authors) > 3:
                author_names += ', et al.'
            results.append({
                'pmid': pmid,
                'title': entry.get('title', ''),
                'authors': author_names,
                'journal': entry.get('source', ''),
                'year': entry.get('pubdate', '')[:4],
                'doi': entry.get('elocationid', ''),
                'url': 'https://pubmed.ncbi.nlm.nih.gov/' + pmid + '/',
            })
        return results
    except Exception as e:
        biosim_log('error', 'pubmed: ' + str(e))
        return None  # None signals error vs empty

def format_citation_bibtex(article):
    """Format article as BibTeX entry."""
    pmid = article.get('pmid', 'unknown')
    first_author = article.get('authors', '').split(',')[0].strip().split(' ')[0]
    year = article.get('year', 'unknown')
    key = first_author + year + pmid[:4]
    return ('@article{' + key + ',\n'
            '  title = {' + article.get('title', '').replace('{', '').replace('}', '') + '},\n'
            '  author = {' + article.get('authors', '') + '},\n'
            '  journal = {' + article.get('journal', '') + '},\n'
            '  year = {' + year + '},\n'
            '  pmid = {' + pmid + '},\n'
            '  url = {' + article.get('url', '') + '},\n'
            '}\n')

def format_citation_ris(article):
    """Format article as RIS (EndNote/Mendeley) format."""
    out = ['TY  - JOUR']
    out.append('PMID  - ' + article.get('pmid', ''))
    out.append('TI  - ' + article.get('title', ''))
    for author in article.get('authors', '').split(','):
        a = author.strip()
        if a:
            out.append('AU  - ' + a)
    out.append('JO  - ' + article.get('journal', ''))
    out.append('PY  - ' + article.get('year', ''))
    out.append('UR  - ' + article.get('url', ''))
    out.append('ER  - ')
    return '\n'.join(out) + '\n'


# ═══════════════════════════════════════════════════
# SYNBIO RESEARCH FUNCTIONS — for scientists
# BioModels import, Hill fitting, stability analysis,
# conservation laws, model reduction
# ═══════════════════════════════════════════════════

def fetch_biomodels(query, max_results=20, timeout=15):
    """Search EBI BioModels database via REST API.
    
    Returns list of dicts with model info, or None on error.
    """
    import urllib.request, urllib.parse, json
    try:
        url = ('https://www.ebi.ac.uk/biomodels/search?query=' +
                urllib.parse.quote(query) +
                '&numResults=' + str(max_results) +
                '&format=json')
        with urllib.request.urlopen(url, timeout=timeout) as f:
            data = json.loads(f.read().decode('utf-8'))
        models = data.get('models', [])
        return [{
            'id': m.get('id', '?'),
            'name': m.get('name', '?'),
            'authors': ', '.join(m.get('submitter', '?').split(',')[:2]) if m.get('submitter') else '',
            'format': m.get('format', '?'),
            'description': (m.get('description') or '')[:200],
            'url': 'https://www.ebi.ac.uk/biomodels/' + m.get('id', ''),
        } for m in models]
    except Exception as e:
        biosim_log('error', 'biomodels: ' + str(e))
        return None

def fetch_biomodel_sbml(model_id, timeout=15):
    """Download SBML for a specific BioModel."""
    import urllib.request
    try:
        url = 'https://www.ebi.ac.uk/biomodels/model/download/' + model_id + '.xml'
        with urllib.request.urlopen(url, timeout=timeout) as f:
            return f.read().decode('utf-8')
    except Exception as e:
        biosim_log('error', 'biomodel download: ' + str(e))
        return None

def fit_hill_equation(x_data, y_data, initial_guess=None):
    """Fit Hill equation y = Vmax * x^n / (K^n + x^n) to data.
    
    Returns dict with Vmax, K, n, R², residuals.
    """
    try:
        from scipy.optimize import curve_fit
        import numpy as np
    except ImportError as _exc:
        biosim_log('warn', 'fit_hill_equation: {}'.format(_exc))
        return None
    
    def hill(x, Vmax, K, n):
        return Vmax * (x ** n) / ((K ** n) + (x ** n))
    
    x_arr = np.array(x_data, dtype=float)
    y_arr = np.array(y_data, dtype=float)
    
    if initial_guess is None:
        initial_guess = [max(y_arr), np.median(x_arr), 1.0]
    
    try:
        popt, pcov = curve_fit(hill, x_arr, y_arr, p0=initial_guess,
                                 bounds=([0, 0.001, 0.1], [1e6, 1e6, 10]),
                                 maxfev=5000)
        Vmax, K, n = popt
        # R²
        y_pred = hill(x_arr, *popt)
        ss_res = np.sum((y_arr - y_pred) ** 2)
        ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        # Standard errors
        try:
            perr = np.sqrt(np.diag(pcov))
        except Exception as _exc:
            biosim_log('warn', 'fit_hill_equation: {}'.format(_exc))
            perr = [0, 0, 0]
        return {
            'Vmax': Vmax, 'K': K, 'n': n,
            'Vmax_se': perr[0], 'K_se': perr[1], 'n_se': perr[2],
            'R2': r2,
            'y_pred': y_pred.tolist(),
        }
    except Exception as e:
        biosim_log('error', 'hill fit: ' + str(e))
        return None

def fit_michaelis_menten(x_data, y_data):
    """Fit Michaelis-Menten: v = Vmax*S / (Km + S)."""
    try:
        from scipy.optimize import curve_fit
        import numpy as np
    except ImportError as _exc:
        biosim_log('warn', 'fit_michaelis_menten: {}'.format(_exc))
        return None
    
    def mm(s, Vmax, Km):
        return Vmax * s / (Km + s)
    
    x_arr = np.array(x_data, dtype=float)
    y_arr = np.array(y_data, dtype=float)
    
    try:
        popt, pcov = curve_fit(mm, x_arr, y_arr,
                                 p0=[max(y_arr), np.median(x_arr)],
                                 bounds=([0, 0.001], [1e6, 1e6]),
                                 maxfev=5000)
        Vmax, Km = popt
        y_pred = mm(x_arr, *popt)
        ss_res = np.sum((y_arr - y_pred) ** 2)
        ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        try:
            perr = np.sqrt(np.diag(pcov))
        except Exception as _exc:
            biosim_log('warn', 'fit_michaelis_menten: {}'.format(_exc))
            perr = [0, 0]
        return {
            'Vmax': Vmax, 'Km': Km,
            'Vmax_se': perr[0], 'Km_se': perr[1],
            'R2': r2,
            'y_pred': y_pred.tolist(),
        }
    except Exception as e:
        biosim_log('error', 'mm fit: ' + str(e))
        return None

def compute_jacobian_numeric(f, x0, t=0, eps=1e-7):
    """Compute Jacobian matrix numerically at point x0.
    
    f(x, t) returns rate vector dx/dt.
    Returns J[i][j] = ∂f_i/∂x_j
    """
    try:
        import numpy as np
        x0 = np.array(x0, dtype=float)
        n = len(x0)
        J = np.zeros((n, n))
        f0 = np.array(f(x0, t), dtype=float)
        for j in range(n):
            xp = x0.copy()
            xp[j] += eps
            fp = np.array(f(xp, t), dtype=float)
            J[:, j] = (fp - f0) / eps
        return J
    except Exception as _exc:
        biosim_log('warn', 'compute_jacobian_numeric: {}'.format(_exc))
        return None

def eigenvalues_stability(J):
    """Analyze stability from Jacobian eigenvalues.
    
    Returns dict with eigenvalues, max_real, type (stable/unstable/saddle/etc).
    """
    try:
        import numpy as np
        eigs = np.linalg.eigvals(J)
        real_parts = eigs.real
        imag_parts = eigs.imag
        max_re = float(real_parts.max())
        min_re = float(real_parts.min())
        
        has_imag = any(abs(im) > 1e-9 for im in imag_parts)
        
        if max_re < -1e-9:
            stability = 'stable_focus' if has_imag else 'stable_node'
        elif max_re > 1e-9:
            if min_re < -1e-9:
                stability = 'saddle'
            else:
                stability = 'unstable_focus' if has_imag else 'unstable_node'
        else:
            stability = 'marginal'
        
        return {
            'eigenvalues': [(float(e.real), float(e.imag)) for e in eigs],
            'max_real': max_re,
            'min_real': min_re,
            'has_complex': has_imag,
            'stability': stability,
            'is_stable': max_re < 0,
        }
    except Exception as e:
        biosim_log('error', 'stability: ' + str(e))
        return None

def find_conservation_laws(stoichiometry_matrix):
    """Find conservation laws (left null space of stoichiometry matrix).
    
    A conserved moiety c satisfies c^T * S = 0
    Returns list of conservation vectors.
    """
    try:
        import numpy as np
        S = np.array(stoichiometry_matrix, dtype=float)
        # SVD-based null space
        u, sv, vt = np.linalg.svd(S.T)
        tol = max(S.shape) * sv[0] * 1e-10 if len(sv) > 0 else 1e-10
        rank = (sv > tol).sum() if len(sv) > 0 else 0
        # Null space columns
        null_space = u[:, rank:]
        conservation_laws = []
        for i in range(null_space.shape[1]):
            v = null_space[:, i]
            # Normalize
            max_abs = abs(v).max()
            if max_abs > 1e-9:
                v = v / max_abs
                conservation_laws.append(v.tolist())
        return conservation_laws
    except Exception as e:
        biosim_log('error', 'conservation: ' + str(e))
        return []


# ═══════════════════════════════════════════════════
# MOLBIO RESEARCH FUNCTIONS
# Multiple alignment, phylogenetics, CRISPR guides,
# plasmid maps, advanced sequence stats
# ═══════════════════════════════════════════════════




class MainWelcomeScreen(tk.Toplevel):
    """Welcome dialog for main VIREON 3D Body Simulator window.
    Shows on first launch (unless dismissed).
    Provides 4 entry-point cards.
    """

    def __init__(self, parent, lang='ru'):
        super().__init__(parent)
        self.parent_main = parent
        self.lang = lang
        self.title('VIREON — ' + ('Добро пожаловать' if lang == 'ru'
                                        else 'Welcome'))
        self.geometry('900x640')
        self.configure(bg=Theme.c('bg'))
        self.transient(parent)
        self.resizable(False, False)

        # Center on screen
        self.update_idletasks()
        w = self.winfo_width()
        h = self.winfo_height()
        x = (self.winfo_screenwidth() // 2) - (w // 2)
        y = (self.winfo_screenheight() // 2) - (h // 2)
        self.geometry('{}x{}+{}+{}'.format(w, h, x, y))

        # ─── Header ───
        header = tk.Frame(self, bg=Theme.c('bg'))
        header.pack(fill='x', pady=(20, 4))
        tk.Label(header, text='VIREON',
                  bg=Theme.c('bg'), fg='#0969da',
                  font=('Arial', 26, 'bold')).pack()
        subtitle = ('Симулятор человеческого тела · Болезни · Лекарства · MolBio'
                     if lang == 'ru' else
                     'Cellular & Synthetic Biology Suite · Diseases · Drugs · MolBio')
        tk.Label(header, text=subtitle,
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 11, 'italic')).pack()

        # ─── 4 entry cards ───
        cards_grid = tk.Frame(self, bg=Theme.c('bg'))
        cards_grid.pack(fill='both', expand=True, padx=30, pady=20)

        cards = [
            {
                'icon': '',
                'color': '#1a7f37',
                'title_ru': 'Изучить болезнь',
                'title_en': 'Explore a disease',
                'desc_ru': 'Выбери из 17 болезней, посмотри 3D тело, симптомы, поражённые органы',
                'desc_en': 'Choose from 17 diseases, see 3D body, symptoms, affected organs',
                'action': self._action_choose_disease,
            },
            {
                'icon': '',
                'color': '#8250df',
                'title_ru': 'Подбор лечения',
                'title_en': 'Select treatment',
                'desc_ru': 'Загрузи готовый препарат или свою SMILES молекулу, симулируй эффект',
                'desc_en': 'Load existing drug or your SMILES molecule, simulate effect',
                'action': self._action_treatment,
            },
            {
                'icon': '',
                'color': '#0969da',
                'title_ru': 'Молекулярная биология',
                'title_en': 'Molecular Biology',
                'desc_ru': 'Sequence Editor, Primers, Restriction, BLAST, Codon optimization',
                'desc_en': 'Sequence Editor, Primers, Restriction, BLAST, Codon optimization',
                'action': self._action_molbio,
            },
            {
                'icon': '',
                'color': '#bc4c00',
                'title_ru': 'Synthetic Biology',
                'title_en': 'Synthetic Biology',
                'desc_ru': 'Конструктор генных схем, ODE/Gillespie, FBA, AI Copilot, COMBINE',
                'desc_en': 'Gene circuit constructor, ODE/Gillespie, FBA, AI Copilot, COMBINE',
                'action': self._action_synbio,
            },
        ]

        for i, card in enumerate(cards):
            row = i // 2
            col = i % 2
            cards_grid.columnconfigure(col, weight=1, minsize=380)
            cards_grid.rowconfigure(row, weight=1)

            cf = tk.Frame(cards_grid, bg=Theme.c('panel'),
                            highlightthickness=2,
                            highlightbackground=Theme.c('border'),
                            cursor='hand2')
            cf.grid(row=row, column=col, padx=10, pady=10, sticky='nsew')

            # Title row with icon
            head = tk.Frame(cf, bg=Theme.c('panel'))
            head.pack(fill='x', padx=18, pady=(16, 4))
            tk.Label(head, text=card['icon'],
                      bg=Theme.c('panel'), fg=card['color'],
                      font=('Arial', 32)).pack(side='left', padx=(0, 8))
            tk.Label(head, text=card['title_' + lang],
                      bg=Theme.c('panel'), fg=card['color'],
                      font=('Arial', 15, 'bold')).pack(side='left', pady=12)

            tk.Label(cf, text=card['desc_' + lang],
                      bg=Theme.c('panel'), fg=Theme.c('text'),
                      font=('Arial', 10), wraplength=320,
                      justify='left').pack(padx=18, pady=(0, 8), anchor='w')

            def make_handler(c):
                def handler(e=None):
                    self.destroy()
                    c['action']()
                return handler

            # Whole card is clickable
            for w in [cf, head] + list(head.winfo_children()):
                w.bind('<Button-1>', make_handler(card))
                w.config(cursor='hand2')

            # Hover effects
            def make_hover(frame, normal='#d0d7de', hover='#0969da'):
                def enter(e):
                    try: frame.config(highlightbackground=hover)
                    except Exception as _exc:
                        biosim_log('warn', 'enter: {}'.format(_exc))
                        pass
                def leave(e):
                    try: frame.config(highlightbackground=normal)
                    except Exception as _exc:
                        biosim_log('warn', 'leave: {}'.format(_exc))
                        pass
                return enter, leave
            enter, leave = make_hover(cf)
            cf.bind('<Enter>', enter)
            cf.bind('<Leave>', leave)

        # ─── Bottom ───
        bottom = tk.Frame(self, bg=Theme.c('bg'))
        bottom.pack(fill='x', pady=14, padx=20)
        self.dont_show_var = tk.IntVar(value=0)
        tk.Checkbutton(bottom,
                        text=('Не показывать при запуске' if lang == 'ru'
                              else 'Don\'t show on startup'),
                        variable=self.dont_show_var,
                        bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                        selectcolor='#ffffff',
                        activebackground='#f5f7fa',
                        font=('Arial', 9)).pack(side='left')

        tk.Button(bottom,
                   text=('Закрыть и продолжить →' if lang == 'ru'
                          else 'Close and continue →'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=6,
                   command=self._close_normally).pack(side='right')

    def _close_normally(self):
        try:
            if self.dont_show_var.get():
                # Save preference
                import json
                settings_path = os.path.expanduser('~/.biosim/settings.json')
                os.makedirs(os.path.dirname(settings_path), exist_ok=True)
                settings = {}
                if os.path.exists(settings_path):
                    try:
                        with open(settings_path) as f:
                            settings = json.load(f)
                    except Exception as _exc:
                        biosim_log('warn', '_close_normally: {}'.format(_exc))
                settings['main_skip_welcome'] = True
                with open(settings_path, 'w') as f:
                    json.dump(settings, f, indent=2)
        except Exception as _exc:
            biosim_log('warn', '_close_normally: {}'.format(_exc))
        self.destroy()

    def _action_choose_disease(self):
        """Focus on disease search."""
        self.destroy()
        try:
            # Highlight search box (find via parent_main)
            if hasattr(self.parent_main, '_search_entry'):
                self.parent_main._search_entry.focus_set()
        except Exception as _exc:
            biosim_log('warn', '_action_choose_disease: {}'.format(_exc))

    def _action_treatment(self):
        self.destroy()

    def _action_molbio(self):
        self.destroy()
        try:
            self.parent_main._open_molbio()
        except Exception as _exc:
            biosim_log('warn', '_action_molbio: {}'.format(_exc))

    def _action_synbio(self):
        self.destroy()
        try:
            self.parent_main._open_synbio()
        except Exception as _exc:
            biosim_log('warn', '_action_synbio: {}'.format(_exc))


class Theme:
    """Central design system for VIREON — colors, fonts, and widget factories.

    Solves the macOS tk.Button bug (where bg is ignored until hover) by always
    setting highlightbackground, and provides a consistent professional look.
    Switch THEME between 'dark' and 'light' (medical white) in one place.
    """
    # ── Active palette (medical light by default; professional, clean) ──
    MODE = 'light'

    PALETTES = {
        # Dark theme — "Graphite Neon" (brightened for better night visibility)
        'dark': {
            'bg':        '#1a1a1f',  # graphite window background
            'panel':     '#26262e',  # panels / cards (brighter than before)
            'panel2':    '#32323c',  # raised elements
            'border':    '#3d3d47',
            'text':      '#f7f7f9',  # primary text (brighter near-white)
            'text_dim':  '#a0a0ad',  # secondary text (brighter muted)
            'accent':    '#2ee0fa',  # neon cyan (brighter)
            'accent_fg': '#101014',
            'success':   '#b4f04a',  # neon lime (brighter)
            'warn':      '#ffc022',
            'danger':    '#ff6b6b',  # brighter red
            'canvas':    '#202026',
        },
        # Light theme — "Ocean" (fresh turquoise, calm scientific look)
        'light': {
            'bg':        '#eaf4f7',  # ocean mist background
            'panel':     '#ffffff',  # clean white cards
            'panel2':    '#dcecf0',  # subtle raised (turquoise tint)
            'border':    '#cce6ec',
            'text':      '#083344',  # deep ocean text
            'text_dim':  '#55707a',  # muted teal-grey
            'accent':    '#0891b2',  # ocean cyan
            'accent_fg': '#ffffff',
            'success':   '#0d9488',  # ocean teal
            'warn':      '#9a6700',
            'danger':    '#cf222e',
            'canvas':    '#ffffff',  # pure white for body image
        },
    }

    # Font family: SF Pro / Helvetica on Mac look professional
    FONT = 'Helvetica'
    SIZES = {'h1': 18, 'h2': 14, 'h3': 12, 'body': 11, 'small': 10, 'tiny': 9}

    @classmethod
    def c(cls, key):
        """Get a color from the active palette."""
        return cls.PALETTES[cls.MODE].get(key, '#000000')

    @classmethod
    def f(cls, size='body', bold=False):
        """Get a font tuple."""
        weight = 'bold' if bold else 'normal'
        return (cls.FONT, cls.SIZES.get(size, 11), weight)

    @classmethod
    def btn(cls, parent, text, command=None, primary=False, **kw):
        """Create a button as a tk.Label (works reliably on macOS, unlike tk.Button).

        Monochrome professional style: white/grey by default, green only if primary.
        Text is ALWAYS visible (no hover-to-reveal bug).
        """
        if primary:
            bg, fg = cls.c('success'), '#ffffff'        # green for key actions
            hover = cls._lighten(cls.c('success'))
        else:
            bg, fg = '#ffffff', '#1c2128'               # white with black text
            hover = '#eef1f5'                            # subtle grey on hover
        lbl = tk.Label(parent, text=text, bg=bg, fg=fg,
                       font=kw.pop('font', cls.f('small', bold=True)),
                       cursor='hand2', padx=kw.pop('padx', 12), pady=kw.pop('pady', 6),
                       bd=1, relief='solid', highlightthickness=0, **kw)
        # Black border via a thin frame look — use relief solid + dark border color
        lbl.configure(highlightbackground=Theme.c('border'))
        def on_enter(e): lbl.configure(bg=hover)
        def on_leave(e): lbl.configure(bg=bg)
        def on_click(e):
            if command:
                command()
        lbl.bind('<Enter>', on_enter)
        lbl.bind('<Leave>', on_leave)
        lbl.bind('<Button-1>', on_click)
        return lbl

    @classmethod
    def button(cls, parent, text, command=None, kind='default', **kw):
        """Create a properly-styled button that works on macOS.

        kind: 'default' | 'primary' | 'success' | 'danger' | 'ghost'
        """
        styles = {
            'default': (cls.c('panel2'), cls.c('text')),
            'primary': (cls.c('accent'), cls.c('accent_fg')),
            'success': (cls.c('success'), '#ffffff'),
            'danger':  (cls.c('danger'), '#ffffff'),
            'ghost':   (cls.c('panel'), cls.c('text_dim')),
        }
        bg, fg = styles.get(kind, styles['default'])
        opts = dict(
            text=text, command=command,
            bg=bg, fg=fg,
            activebackground=bg, activeforeground=fg,
            highlightbackground=bg,        # ← critical for macOS
            highlightcolor=bg,
            relief='flat', bd=0,
            font=cls.f('small', bold=True),
            cursor='hand2', padx=12, pady=6,
        )
        opts.update(kw)
        btn = tk.Button(parent, **opts)
        # Subtle hover effect
        hover_bg = cls._lighten(bg)
        btn.bind('<Enter>', lambda e: btn.config(bg=hover_bg, highlightbackground=hover_bg))
        btn.bind('<Leave>', lambda e: btn.config(bg=bg, highlightbackground=bg))
        return btn

    @staticmethod
    def _lighten(hex_color, amt=18):
        """Lighten a hex color for hover."""
        try:
            h = hex_color.lstrip('#')
            r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
            r = min(255, r + amt); g = min(255, g + amt); b = min(255, b + amt)
            return '#{:02x}{:02x}{:02x}'.format(r, g, b)
        except Exception as _exc:
            biosim_log('warn', '_lighten: {}'.format(_exc))
            return hex_color

    @classmethod
    def setup_ttk(cls, root):
        """Configure ttk styles for a consistent professional look."""
        style = ttk.Style(root)
        try:
            style.theme_use('clam')  # most customizable base
        except Exception as _exc:
            biosim_log('warn', 'setup_ttk: {}'.format(_exc))
        style.configure('.', background=cls.c('bg'), foreground=cls.c('text'))
        style.configure('TNotebook', background=cls.c('bg'), borderwidth=0)
        style.configure('TNotebook.Tab',
                        background=cls.c('panel2'), foreground=cls.c('text'),
                        padding=(14, 7), font=cls.f('small', bold=True))
        style.map('TNotebook.Tab',
                  background=[('selected', cls.c('panel'))],
                  foreground=[('selected', cls.c('accent'))])
        style.configure('TFrame', background=cls.c('bg'))
        style.configure('TLabel', background=cls.c('bg'), foreground=cls.c('text'))
        return style


class MainApp(tk.Tk):
    def __init__(self, lang='ru'):
        super().__init__()
        self.lang = lang
        self.L = LANG[lang]
        self.cur_dis = None
        self.cur_dat = None
        self.drug_dose = tk.IntVar(value=50)
        self.mut_gene = tk.StringVar()
        self.mut_val = tk.IntVar(value=0)
        self.drug_var = tk.StringVar()
        self.days_var = tk.IntVar(value=100)
        self.dose_mg_var = tk.IntVar(value=100)
        self.weight_var = tk.IntVar(value=70)
        self.height_var = tk.IntVar(value=170)
        self._list_items = []
        self._api_result = None
        # User-added (non-persistent) diseases for this session
        self._custom_diseases = {}

        # Load saved theme preference (light by default) BEFORE building UI
        self._load_theme_preference()

        self.title(self.L['title'])
        self.geometry('1320x900')
        Theme.setup_ttk(self)
        self.configure(bg=Theme.c('bg'))
        self.minsize(1000, 700)
        self._build()

        # Ensure the whole program quits cleanly when the window is closed
        # (otherwise child windows / threads can keep the process alive).
        self.protocol('WM_DELETE_WINDOW', self._on_quit)

        # Show welcome screen if not dismissed
        self.after(300, self._maybe_show_welcome)
        # Apply the modern flat button/font style across the WHOLE app, after
        # all tabs are built. Re-run periodically-ish on tab switches via the
        # public helper below.
        self.after(400, self._global_restyle)

    def _global_restyle(self):
        """Walk every widget in the app and apply the flat modern style
        (Ocean / Graphite-Neon): remove grey Windows borders from buttons,
        switch Arial → modern UI font on labels/buttons. Safe to call repeatedly.
        """
        accent_like = {'#1a7f37', '#0b7a4b', '#0d9488', '#0550ae', '#0969da',
                       '#2f81f7', '#0891b2', '#8250df', '#6639ba'}
        danger_like = {'#cf222e', '#d1242f', '#a40e26'}
        try:
            font_main = BIOSIM_FONT_FALLBACK
        except Exception:
            font_main = 'Helvetica'
        panel_bg = biosim_theme('bg')
        def walk(w):
            for c in w.winfo_children():
                cls = c.__class__.__name__
                try:
                    if cls == 'Button':
                        cur = str(c.cget('bg')).lower()
                        if cur in danger_like:
                            style_button(c, 'danger')
                        elif cur in accent_like:
                            style_button(c, 'primary')
                        else:
                            style_button(c, 'normal')
                        # GUARANTEE readable text: pick dark/light by the button's
                        # actual background luminance (fixes white-on-white and
                        # white-on-light in both themes).
                        try:
                            bgc = str(c.cget('bg'))
                            c.configure(fg=readable_text_on(bgc))
                        except Exception:
                            pass
                        try:
                            if 'Arial' in str(c.cget('font')):
                                c.configure(font=(font_main, 10))
                        except Exception:
                            pass
                    elif cls == 'Label':
                        try:
                            # Fix invisible label text: if the label's text colour
                            # is (near) the same brightness as its background, recolour.
                            lbg = str(c.cget('bg')) or panel_bg
                            lfg = str(c.cget('fg'))
                            # Only correct near-white text on light bg / near-dark on dark.
                            if lfg.lower() in ('#ffffff', 'white', '#f5f7fa', '#f5f5f5',
                                               '#e6edf3', '#f4f4f5'):
                                import re as _re2
                                if _hex_luminance(lbg if lbg.startswith('#') else panel_bg) > 0.55:
                                    c.configure(fg=biosim_theme('text'))
                            fnt = str(c.cget('font'))
                            if 'Arial' in fnt:
                                size = 10; weight = 'normal'
                                import re as _re
                                m = _re.search(r"(\d+)", fnt)
                                if m: size = int(m.group(1))
                                if 'bold' in fnt: weight = 'bold'
                                c.configure(font=(font_main, size, weight))
                        except Exception:
                            pass
                except Exception:
                    pass
                walk(c)
        try:
            walk(self)
        except Exception:
            pass

    def _on_quit(self):
        """Cleanly shut down the entire application."""
        try:
            # Close all child windows (Toplevels) first
            for w in self.winfo_children():
                if isinstance(w, tk.Toplevel):
                    try:
                        w.destroy()
                    except Exception as _exc:
                        biosim_log('warn', '_on_quit: {}'.format(_exc))
        except Exception as _exc:
            biosim_log('warn', '_on_quit: {}'.format(_exc))
        try:
            self.quit()      # stop the mainloop
        except Exception as _exc:
            biosim_log('warn', '_on_quit: {}'.format(_exc))
        try:
            self.destroy()   # destroy the window
        except Exception as _exc:
            biosim_log('warn', '_on_quit: {}'.format(_exc))
        # Force the process to exit so no Python stays running in the background
        import os
        os._exit(0)


    def _load_theme_preference(self):
        """Read the saved theme (light/dark) from settings and apply it to
        the Theme class before the UI is built."""
        try:
            import json
            settings_path = os.path.expanduser('~/.biosim/settings.json')
            if os.path.exists(settings_path):
                with open(settings_path) as f:
                    settings = json.load(f)
                mode = settings.get('theme_mode', 'light')
                if mode in ('light', 'dark'):
                    Theme.MODE = mode
        except Exception as _exc:
            biosim_log('warn', '_load_theme_preference: {}'.format(_exc))
            Theme.MODE = 'light'
        # Sync the shared theme (SynBio/MolBio) with the loaded mode
        try:
            set_biosim_theme(Theme.MODE)
        except Exception as _exc:
            biosim_log('warn', '_load_theme_preference: {}'.format(_exc))

    def _save_theme_preference(self):
        """Persist the current theme choice to settings."""
        try:
            import json
            os.makedirs(os.path.expanduser('~/.biosim'), exist_ok=True)
            settings_path = os.path.expanduser('~/.biosim/settings.json')
            settings = {}
            if os.path.exists(settings_path):
                with open(settings_path) as f:
                    settings = json.load(f)
            settings['theme_mode'] = Theme.MODE
            with open(settings_path, 'w') as f:
                json.dump(settings, f, indent=2)
        except Exception as _exc:
            biosim_log('warn', '_save_theme_preference: {}'.format(_exc))

    def toggle_theme(self):
        """Switch between light and dark theme, save the choice, and rebuild
        the UI so the change is visible immediately. User data is NOT touched —
        only colours change."""
        Theme.MODE = 'dark' if Theme.MODE == 'light' else 'light'
        self._save_theme_preference()
        # Sync the shared theme used by SynBio (part2) and MolBio (part4)
        try:
            set_biosim_theme(Theme.MODE)
        except Exception as _exc:
            biosim_log('warn', 'toggle_theme: {}'.format(_exc))
        # Rebuild the main window UI with the new colours. Open child windows
        # (SynBio/MolBio) keep their own state; they pick up the theme when
        # reopened. We only rebuild the main app body here.
        try:
            self.configure(bg=Theme.c('bg'))
            Theme.setup_ttk(self)
            # Destroy and rebuild the body
            for w in self.winfo_children():
                try:
                    w.destroy()
                except Exception as _exc:
                    biosim_log('warn', 'toggle_theme: {}'.format(_exc))
            self._build()
            # Re-apply the flat modern button style with the NEW theme colours.
            self.after(200, self._global_restyle)
            is_ru = (self.lang == 'ru')
            mode_name = ('тёмная' if Theme.MODE == 'dark' else 'светлая') if is_ru \
                        else ('dark' if Theme.MODE == 'dark' else 'light')
            # Brief non-blocking confirmation
            try:
                self.title(self.L['title'] + '  ·  ' +
                           ('Тема: ' if is_ru else 'Theme: ') + mode_name)
            except Exception as _exc:
                biosim_log('warn', 'toggle_theme: {}'.format(_exc))
        except Exception as e:
            biosim_log('error', 'theme toggle: ' + str(e))

    def _maybe_show_welcome(self):
        """Show MainWelcomeScreen unless user said don't show."""
        try:
            import json
            settings_path = os.path.expanduser('~/.biosim/settings.json')
            if os.path.exists(settings_path):
                with open(settings_path) as f:
                    settings = json.load(f)
                if settings.get('main_skip_welcome'):
                    return
            MainWelcomeScreen(self, self.lang)
        except Exception as e:
            biosim_log('error', 'welcome screen: ' + str(e))

    def _frm(self, parent, title):
        return tk.LabelFrame(parent, text=title, bg=Theme.c('panel'), fg='#0969da',
                              font=('Arial', 11, 'bold'), relief='flat', bd=1,
                              highlightbackground=Theme.c('border'))

    def _build(self):
        # ─── Menu bar (always visible at top) ───
        is_ru = (self.lang == 'ru')
        menubar = tk.Menu(self)
        tools_menu = tk.Menu(menubar, tearoff=0)
        tool_items = [
            ('' + ('Профиль пациента' if is_ru else 'Patient Profile'), self._open_patient_profile),
            ('' + ('ИИ-ассистент' if is_ru else 'AI Assistant'), self._open_ai_assistant),
            ('Lab values', self._open_lab_values),
            ('' + ('Сравнить болезни' if is_ru else 'Compare'), self._open_disease_comparison),
            ('' + ('Взаимодействия' if is_ru else 'Interactions'), self._open_drug_interactions),
            ('' + ('Отчёт' if is_ru else 'Report'), self._open_report_generator),
            ('' + ('Журнал' if is_ru else 'Notes'), self._open_notes),
        ]
        for label, cmd in tool_items:
            tools_menu.add_command(label=label, command=cmd)
        menubar.add_cascade(label='' + ('Инструменты' if is_ru else 'Tools'),
                             menu=tools_menu)

        stats_menu = tk.Menu(menubar, tearoff=0)
        for label, cmd in [
            ('PK/PD', self._open_pk_pd),
            ('' + ('Фармакология PRO (PBPK)' if is_ru else 'Pharmacology PRO (PBPK)'), self._open_pharmacology_pro),
            ('' + ('Фармакокинетика PRO' if is_ru else 'Advanced PK'), self._open_pk_advanced),
            ('' + ('Калькулятор доз' if is_ru else 'Dose Calculator'), self._open_dose_calculator),
            ('Survival (Kaplan-Meier)', self._open_survival_analysis),
            ('' + ('Статистика' if is_ru else 'Statistics'), self._open_statistics),
            ('' + ('Когорта' if is_ru else 'Cohort'), self._open_cohort_simulation),
            ('' + ('Импорт CSV' if is_ru else 'Import CSV'), self._open_data_import),
            ('Manhattan (GWAS)', self._open_manhattan_plot),
        ]:
            stats_menu.add_command(label=label, command=cmd)
        menubar.add_cascade(label='' + ('Анализ' if is_ru else 'Analysis'),
                             menu=stats_menu)

        science_menu = tk.Menu(menubar, tearoff=0)
        for label, cmd in [
            # MolBio PRO moved to the Molecular Biology tab; SynBio tools moved
            # to the Constructor — they belong to those specific workspaces.
            ('' + ('Эпидемиология' if is_ru else 'Epidemiology'), self._open_epidemiology),
            ('' + ('Эпидемиология PRO' if is_ru else 'Epidemiology PRO'), self._open_epidemiology_pro),
            ('' + ('Эпидемиология PRO+ (SEIRS, сети)' if is_ru else 'Epidemiology PRO+ (SEIRS, networks)'), self._open_epidemiology_advanced),
            ('' + ('Нейрон (1 клетка)' if is_ru else 'Neuron (single)'), self._open_neuroscience),
            ('' + ('Сеть нейронов' if is_ru else 'Neural network'), self._open_neuro_network),
            ('' + ('Нейробиология PRO+ (типы, память, обучение)' if is_ru else 'Neuroscience PRO+ (types, memory, learning)'), self._open_neuroscience_advanced),
            ('' + ('Пластичность (STDP)' if is_ru else 'Plasticity (STDP)'), self._open_neuro_plasticity),
            ('' + ('Дерево NJ' if is_ru else 'NJ Tree'), self._open_phylo_builder),
            ('' + ('Молекулярный докинг' if is_ru else 'Molecular Docking'), self._open_docking),
            ('' + ('Выравнивание MSA' if is_ru else 'MSA Alignment'), self._open_msa),
            ('ADME', self._open_adme),
        ]:
            science_menu.add_command(label=label, command=cmd)
        menubar.add_cascade(label='' + ('Наука' if is_ru else 'Science'),
                             menu=science_menu)

        # ── View menu: show/hide panels with checkmarks ──
        view_menu = tk.Menu(menubar, tearoff=0)
        self._panel_vars = {}
        panel_labels = [
            ('left', 'Поиск и информация' if is_ru else 'Search & Info'),
            ('center', 'График и тело' if is_ru else 'Graph & Body'),
            ('right', 'Результаты и управление' if is_ru else 'Results & Controls'),
        ]
        for pkey, plabel in panel_labels:
            var = tk.BooleanVar(value=True)
            self._panel_vars[pkey] = var
            view_menu.add_checkbutton(
                label=plabel, variable=var,
                command=lambda k=pkey: self._toggle_panel(k))
        menubar.add_cascade(label=('' if is_ru else '') + ('Вид' if is_ru else 'View'),
                             menu=view_menu)

        other_menu = tk.Menu(menubar, tearoff=0)
        for label, cmd in [
            ('' + ('Викторина' if is_ru else 'Quiz'), self._open_anatomy_quiz),
            ('' + ('Настройки' if is_ru else 'Settings'), self._open_settings),
            ('ℹ ' + ('О программе' if is_ru else 'About'), self._open_about),
            ('' + ('Все инструменты' if is_ru else 'All tools'), self._open_all_tools),
        ]:
            other_menu.add_command(label=label, command=cmd)
        menubar.add_cascade(label='' + ('Ещё' if is_ru else 'More'),
                             menu=other_menu)
        self.config(menu=menubar)

        top = tk.Frame(self, bg=Theme.c('panel'), height=44)
        top.pack(fill='x'); top.pack_propagate(False)
        tk.Label(top, text='VIREON', font=('Arial', 17, 'bold'),
                 bg=Theme.c('panel'), fg='#0969da').pack(side='left', padx=16, pady=6)
        tk.Label(top, text=self.L['title'], font=('Arial', 10),
                 bg=Theme.c('panel'), fg=Theme.c('text_dim')).pack(side='left', pady=6)

        # Resizable panels via PanedWindow — drag the dividers to resize.
        main = tk.PanedWindow(self, orient='horizontal', bg='#d0d7de',
                              sashwidth=6, sashrelief='flat', bd=0,
                              showhandle=False)
        main.pack(fill='both', expand=True, padx=6, pady=6)
        self._main_container = main

        # Panels are tracked so they can be shown/hidden from the View menu
        self._panels = {}
        self._panel_order = ['left', 'center', 'right']

        left = tk.Frame(main, bg=Theme.c('bg'), width=290)
        self._build_search(left); self._build_info(left)
        self._panels['left'] = left

        center = tk.Frame(main, bg=Theme.c('bg'))
        self._build_graph(center); self._build_3d(center)
        self._panels['center'] = center

        right = tk.Frame(main, bg=Theme.c('bg'), width=360)
        self._build_text(right); self._build_controls(right)
        self._panels['right'] = right

        # Add to paned window (left & right keep size, center stretches)
        main.add(left, minsize=200, width=290, stretch='never')
        main.add(center, minsize=300, stretch='always')
        main.add(right, minsize=250, width=360, stretch='never')

    def _build_search(self, p):
        f = self._frm(p, self.L['search']); f.pack(fill='x', pady=(0, 4))
        row = tk.Frame(f, bg=Theme.c('panel')); row.pack(fill='x', padx=6, pady=6)
        self.sv = tk.StringVar()
        e = tk.Entry(row, textvariable=self.sv, bg=Theme.c('panel2'), fg=Theme.c('text'),
                     insertbackground='#1c2128', relief='flat', font=('Arial', 11),
                     highlightthickness=1, highlightbackground=Theme.c('border'))
        e.pack(side='left', fill='x', expand=True, ipady=5)
        e.bind('<Return>', lambda ev: self._do_search())
        tk.Button(row, text=self.L['find'], bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                  relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                  padx=8, command=self._do_search).pack(side='right', padx=(4, 2))
        tk.Button(row, text=self.L['add_dis'], bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                  relief='solid', bd=1, font=('Arial', 14, 'bold'), cursor='hand2',
                  width=2, command=self._open_new_disease).pack(side='right')

        self.search_lbl = tk.Label(f, text='', bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                                    font=('Arial', 9))
        self.search_lbl.pack(anchor='w', padx=6)

        self.lb = tk.Listbox(f, bg=Theme.c('panel2'), fg=Theme.c('text'),
                              selectbackground='#0969da', font=('Arial', 11),
                              relief='flat', height=8, highlightthickness=0)
        self.lb.pack(fill='x', padx=6, pady=(0, 6))
        self.lb.bind('<<ListboxSelect>>', self._on_sel)
        self._fill_local('')

    def _fill_local(self, q):
        q = q.lower()
        items = []
        for key, d in DB.items():
            if not q or q in d['en'].lower() or q in d['ru'].lower():
                items.append({'key': key, 'en': d['en'], 'ru': d['ru'], 'src': 'local'})
        for key, d in self._custom_diseases.items():
            if not q or q in key.lower():
                items.append({'key': key, 'en': key, 'ru': key, 'src': 'custom'})
        self._list_items = items
        self.lb.delete(0, 'end')
        for it in items:
            label = it['ru'] if self.lang == 'ru' else it['en']
            tag = ' [custom]' if it['src'] == 'custom' else ''
            self.lb.insert('end', label + tag)

    def _do_search(self):
        q = self.sv.get().strip()
        if not q:
            self._fill_local(''); return
        if len(q) < 2:
            self._fill_local(q); return
        self.search_lbl.config(text=self.L['searching'])
        self.lb.delete(0, 'end')
        self.lb.insert('end', self.L['searching'])

        def worker():
            results = search_diseases(q, DB, max_res=80)
            q_low = q.lower()
            for key in self._custom_diseases:
                if q_low in key.lower():
                    results.insert(0, {'key': key, 'en': key, 'ru': key, 'src': 'custom'})
            def update():
                self._list_items = results
                self.lb.delete(0, 'end')
                if not results:
                    self.lb.insert('end', self.L['not_found'])
                    self.search_lbl.config(text=''); return
                for it in results:
                    label = it.get('ru', '') or it.get('en', '')
                    src = it.get('src', '')
                    badge = ' [{}]'.format(src) if src and src != 'local' else ''
                    self.lb.insert('end', label + badge)
                self.search_lbl.config(text='{}: {}'.format(self.L['found'], len(results)))
            self.after(0, update)

        threading.Thread(target=worker, daemon=True).start()

    def _on_sel(self, _):
        sel = self.lb.curselection()
        if not sel:
            return
        idx = int(sel[0])
        if idx >= len(self._list_items):
            return
        item = self._list_items[idx]
        key = item.get('key', '')
        # PubMed article — show in info panel
        if item.get('src') == 'PubMed':
            self._show_info_pubmed(item)
            return
        if key in DB:
            self.cur_dis = key
            self.cur_dat = DB[key]
            self._show_info_local(key)
            self._update_menus(key)
        elif key in self._custom_diseases:
            self.cur_dis = key
            self.cur_dat = self._custom_diseases[key]
            self._show_info_custom(key)
        else:
            self.cur_dis = key
            self.cur_dat = None
            self._show_info_api_item(item)

    def _show_info_pubmed(self, item):
        """Display a PubMed article in the info panel."""
        is_ru = (self.lang == 'ru')
        lines = [
            ('PubMed\n', 'hdr'),
            ('-' * 26 + '\n', 'dim'),
            (item.get('en', '') + '\n\n', 'norm'),
        ]
        if item.get('authors'):
            lines.append((('Авторы: ' if is_ru else 'Authors: ') + item['authors'] + '\n', 'lbl'))
        if item.get('journal'):
            lines.append((('Журнал: ' if is_ru else 'Journal: ') + item['journal'] + '\n', 'norm'))
        if item.get('year'):
            lines.append((('Год: ' if is_ru else 'Year: ') + item['year'] + '\n', 'norm'))
        if item.get('pmid'):
            lines.append(('PMID: ' + item['pmid'] + '\n', 'dim'))
        if item.get('url'):
            lines.append(('\n' + item['url'] + '\n', 'lbl'))
        lines.append(('\n' + ('Двойной клик — открыть в браузере' if is_ru
                                    else 'Double-click to open in browser') + '\n', 'dim'))
        self._iw(lines)
        self._pubmed_url = item.get('url', '')
        try:
            self.info_txt.bind('<Double-Button-1>', self._open_pubmed_link)
        except Exception as _exc:
            biosim_log('warn', '_show_info_pubmed: {}'.format(_exc))

    def _open_pubmed_link(self, event=None):
        url = getattr(self, '_pubmed_url', '')
        if url:
            try:
                import webbrowser
                webbrowser.open(url)
            except Exception as _exc:
                biosim_log('warn', '_open_pubmed_link: {}'.format(_exc))

    def _open_new_disease(self):
        NewDiseaseDialog(self, self.lang, self._handle_new_disease)

    def _open_synbio(self):
        """Open Synthetic Biology window (independent, non-modal)."""
        if hasattr(self, '_synbio_window') and self._synbio_window is not None:
            try:
                if self._synbio_window.winfo_exists():
                    self._synbio_window.lift()
                    self._synbio_window.focus_force()
                    return
            except Exception as _exc:
                biosim_log('warn', '_open_synbio: {}'.format(_exc))
        self._synbio_window = SynBioWindow(self, self.lang)

    def _open_molbio(self):
        """Open Molecular Biology workspace (independent, non-modal)."""
        if hasattr(self, '_molbio_window') and self._molbio_window is not None:
            try:
                if self._molbio_window.winfo_exists():
                    self._molbio_window.lift()
                    self._molbio_window.focus_force()
                    self._molbio_load_pending()
                    return
            except Exception as _exc:
                biosim_log('warn', '_open_molbio: {}'.format(_exc))
        try:
            self._molbio_window = MolBioWindow(self, self.lang)
            # If a disease (or other section) left a sequence waiting, load it now
            # so the window isn't empty.
            self.after(300, self._molbio_load_pending)
        except Exception as e:
            biosim_log('error', 'MolBio failed to open: ' + str(e))
            messagebox.showerror('VIREON', 'MolBio error: ' + str(e))

    def _molbio_load_pending(self):
        """Load any sequence that was queued for MolBio (e.g. from a disease
        applied across all workspaces) into the open MolBio window."""
        seq = getattr(self, '_molbio_pending_sequence', '')
        if not seq:
            return
        mw = getattr(self, '_molbio_window', None)
        if mw is None:
            return
        try:
            if mw.winfo_exists() and hasattr(mw, 'load_external_sequence'):
                mw.load_external_sequence(seq)
                # clear so it doesn't reload every time
                self._molbio_pending_sequence = ''
        except Exception as _exc:
            biosim_log('warn', '_molbio_load_pending: {}'.format(_exc))

    def _handle_new_disease(self, name, desc, genes_text, symptoms_text, organs_text, dis_type):
        self._rw([(self.L['loading'] + '...\n', 'wrn')])
        self.update()
        def worker():
            def cb(s, t, m):
                self.after(0, lambda: self._rw(
                    [(self.L['loading'] + ':\n', 'wrn'),
                     ('[{}/{}] {}\n'.format(s, t, m), 'dim')]))
            R = analyze_custom_disease(name, desc, genes_text, symptoms_text, organs_text, dis_type, cb)
            genes_list = [g.strip().upper() for g in genes_text.split(',') if g.strip()]
            organs_list = [o.strip() for o in organs_text.split(',') if o.strip()]
            self._custom_diseases[name] = {
                'en': name, 'ru': name,
                'genes': genes_list,
                'proteins': [p.get('name', '') for p in R.get('proteins', [])],
                'organs_ru': organs_list, 'organs_en': organs_list,
                'organ_map': organs_list,
                'desc_ru': desc, 'desc_en': desc,
                'drugs': {},
                'params': {'base': 0.7, 'prot': 0.3, 'rec': 0.05},
                'surv': 60, 'recov': 30, 'mgmt': 50,
                'type': dis_type, 'symptoms': R.get('symptoms', []),
                'methods': R.get('methods', []),
                'pathways': R.get('pathways', []),
                'publications': R.get('publications', []),
            }
            self.after(0, lambda: self._show_custom_disease(name))
            self.after(0, lambda: self._fill_local(''))
        threading.Thread(target=worker, daemon=True).start()

    def _show_custom_disease(self, name):
        d = self._custom_diseases[name]
        self.cur_dis = name
        self.cur_dat = d
        self._show_info_custom(name)
        self._show_methods(d.get('methods', []))

    def _show_info_custom(self, name):
        d = self._custom_diseases[name]
        lines = [
            (name + '\n', 'hdr'),
            ('-' * 26 + '\n', 'dim'),
            ('Type: ' + d.get('type', '') + '\n', 'dim'),
            (d.get('desc_en', '') + '\n\n', 'norm'),
        ]
        if d.get('genes'):
            lines.append((self.L['genes'] + ':\n', 'lbl'))
            for g in d['genes']:
                lines.append(('  ' + g + '\n', 'gene'))
        if d.get('symptoms'):
            lines.append(('Symptoms:\n', 'lbl'))
            for s in d['symptoms']:
                lines.append(('  ' + s + '\n', 'norm'))
        if d.get('organs_en'):
            lines.append((self.L['organs'] + ':\n', 'lbl'))
            for o in d['organs_en']:
                lines.append(('  ' + o + '\n', 'org'))
        self._iw(lines)

    def _show_methods(self, methods):
        if not methods:
            self._rw([(self.L['not_found'] + '\n', 'dim')])
            return
        lines = [(self.L['candidates'] + '\n', 'ttl'),
                  ('-' * 26 + '\n', 'dim'),
                  ('Use one of the molecules below:\n\n', 'dim')]
        for i, m in enumerate(methods, 1):
            lines.append(('METHOD {}: {}\n'.format(i, m.get('type', 'Inhibitor')), 'ttl'))
            lines.append(('  Target: {}\n'.format(m['target'][:50]), 'nrm'))
            lines.append(('  Gene:   {}\n'.format(m['gene']), 'gene'))
            lines.append(('  ' + m.get('target_function', '')[:120] + '\n', 'dim'))
            lines.append(('  ' + self.L['molecules'] + ':\n', 'lbl'))
            for dr in m.get('drugs', [])[:5]:
                ic = dr.get('ic50_nm')
                ic_s = '  IC50={} nM'.format(ic) if ic else ''
                tag = 'ok' if dr['pchembl'] > 7 else 'wrn'
                lines.append(('    - {} (pChEMBL={:.1f}{})\n'.format(
                    dr['mol_name'], dr['pchembl'], ic_s), tag))
            lines.append(('\n', 'nrm'))
        lines.append(('To use: type drug name in "Drug analysis" field below.\n', 'dim'))
        self._rw(lines)


    def _show_info_api_item(self, item):
        """Handle disease selected from API search (not in local DB).
        Shows description immediately, then fetches genes via UniProt in background."""
        is_ru = (self.lang == 'ru')
        name = item.get('ru', '') or item.get('en', '') or item.get('key', '')
        desc = item.get('desc', '')
        src = item.get('src', '')

        # Detect probable organs from description
        organs_detected = []
        desc_low = (desc + ' ' + name).lower()
        organ_keywords = {
            'Lungs':    ['lung','pulmon','respirat','лёг','лег','дыха','пневмон'],
            'Heart':    ['heart','cardi','коронар','серд','миокард'],
            'Liver':    ['liver','hepati','цирроз','печен','гепат'],
            'Kidneys':  ['kidney','renal','нефро','почек','почк','ренал'],
            'Brain':    ['brain','cerebr','neuro','мозг','нерв','нейр','энцефал'],
            'Pancreas': ['pancrea','diabet','подж','диабет'],
            'Muscles':  ['muscle','skin','derma','melanoma','кожи','кож','меланом','мышц','грудн'],
            'Blood':    ['blood','hemo','anemia','кров','гемор','анеми','вирус','virus'],
        }
        for organ, keywords in organ_keywords.items():
            if any(kw in desc_low for kw in keywords):
                organs_detected.append(organ)

        # Build a usable disease record IMMEDIATELY (works even if gene fetch fails)
        disease_record = {
            'en': name, 'ru': name,
            'genes': list(item.get('genes', [])),
            'proteins': [],
            'organs_ru': organs_detected, 'organs_en': organs_detected,
            'organ_map': organs_detected,
            'desc_ru': desc, 'desc_en': desc,
            'drugs': {},
            'params': {'base': 0.7, 'prot': 0.3, 'rec': 0.05},
            'surv': 60, 'recov': 30, 'mgmt': 50,
            'src': src,
        }
        self._custom_diseases[name] = disease_record
        self.cur_dis = name
        self.cur_dat = disease_record

        # Show info right away (don't wait for network)
        def render(genes, loading=False):
            lines = [(name + '\n', 'hdr'), ('-' * 26 + '\n', 'dim'),
                     ((('Источник: ' if is_ru else 'Source: ') + src + '\n\n'), 'dim')]
            if desc:
                lines.append((desc + '\n\n', 'norm'))
            else:
                lines.append((('Описание недоступно из источника.\n\n' if is_ru
                               else 'No description available from source.\n\n'), 'dim'))
            if genes:
                lines.append((self.L['genes'] + ' (UniProt):\n', 'lbl'))
                for g in genes:
                    lines.append(('  ' + g + '\n', 'gene'))
                lines.append(('\n', 'norm'))
            elif loading:
                lines.append((('Загрузка генов из UniProt...\n\n' if is_ru
                               else 'Loading genes from UniProt...\n\n'), 'wrn'))
            else:
                lines.append((('Специфичные гены не найдены — будут использованы шаблоны.\n\n'
                               if is_ru else
                               'No specific genes found — drug templates will be used.\n\n'), 'dim'))
            if organs_detected:
                lines.append((self.L['organs'] + (' (определены)' if is_ru else ' (detected)') + ':\n', 'lbl'))
                for o in organs_detected:
                    lines.append(('  ' + o + '\n', 'org'))
                lines.append(('\n', 'norm'))
            lines.append((('Готово! Теперь можно:\n' if is_ru else 'Ready! You can now:\n'), 'ttl'))
            lines.append((('  • Создать лекарство\n' if is_ru else '  - Create a drug\n'), 'norm'))
            lines.append((('  • Запустить симуляцию\n' if is_ru else '  - Run simulation\n'), 'norm'))
            lines.append((('  • Анализировать препараты\n' if is_ru else '  - Analyze drugs\n'), 'norm'))
            self._iw(lines)

        # Initial render (with any genes from search, or loading state)
        render(disease_record['genes'], loading=not disease_record['genes'])

        # Fetch genes in background (with error protection)
        def worker():
            genes = []
            try:
                genes = fetch_disease_genes(name) or []
            except Exception as _exc:
                biosim_log('warn', 'worker: {}'.format(_exc))
                genes = []
            def update():
                if genes:
                    disease_record['genes'] = genes
                    self._custom_diseases[name] = disease_record
                    # Update mutation menu
                    if hasattr(self, 'mut_menu'):
                        dg = list(genes) + [g for g in HEALTHY_GENES if g not in genes]
                        self.mut_menu['values'] = dg
                        if dg:
                            self.mut_menu.current(0)
                            self.mut_gene.set(dg[0])
                # Re-render with final gene info (only if this disease is still selected)
                if self.cur_dis == name:
                    render(genes, loading=False)
            self.after(0, update)

        threading.Thread(target=worker, daemon=True).start()

    def _build_info(self, p):
        f = self._frm(p, self.L['dis_info']); f.pack(fill='both', expand=True)
        self.info_txt = scrolledtext.ScrolledText(
            f, bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 10), relief='flat',
            wrap='word', state='disabled', highlightthickness=0)
        self.info_txt.pack(fill='both', expand=True, padx=4, pady=4)
        for tag, fg, fnt in [
            ('hdr', '#0969da', ('Arial', 11, 'bold')),
            ('h1', '#0969da', ('Arial', 13, 'bold')),
            ('ttl', '#0969da', ('Arial', 11, 'bold')),
            ('lbl', '#57606a', ('Arial', 10, 'bold')),
            ('label', '#1a7f37', ('Arial', 10, 'bold')),
            ('gene', '#1a7f37', ('Arial', 10)),
            ('prot', '#bc4c00', ('Arial', 10)),
            ('org', '#cf222e', ('Arial', 10)),
            ('norm', '#1c2128', ('Arial', 10)),
            ('text', '#1c2128', ('Arial', 10)),
            ('dim', '#57606a', ('Arial', 9)),
            ('wrn', '#9a6700', ('Arial', 10)),
            ('link', '#0550ae', ('Arial', 10, 'underline')),
        ]:
            self.info_txt.tag_config(tag, foreground=fg, font=fnt)
        self._iw([(self.L['hint'] + '\n', 'norm')])

    def _iw(self, lines):
        self.info_txt.config(state='normal')
        self.info_txt.delete('1.0', 'end')
        for txt, tag in lines:
            self.info_txt.insert('end', txt, tag)
        self.info_txt.config(state='disabled')

    def _show_info_local(self, key):
        d = DB[key]
        name = d['ru'] if self.lang == 'ru' else d['en']
        desc = d['desc_ru'] if self.lang == 'ru' else d['desc_en']
        orgs = d['organs_ru'] if self.lang == 'ru' else d['organs_en']
        lines = [(name + '\n', 'hdr'), ('-' * 26 + '\n', 'dim'),
                 (desc + '\n\n', 'norm'), (self.L['genes'] + ':\n', 'lbl')]
        for g in d['genes']:
            lines.append(('  ' + g + '\n', 'gene'))
        lines.append((self.L['proteins'] + ':\n', 'lbl'))
        for pr in d['proteins']:
            lines.append(('  ' + pr + '\n', 'prot'))
        lines.append((self.L['organs'] + ':\n', 'lbl'))
        for o in orgs:
            lines.append(('  ' + o + '\n', 'org'))
        self._iw(lines)

    def _build_graph(self, p):
        f = self._frm(p, self.L['graph'])
        f.pack(fill='both', expand=True, pady=(0, 4))
        self.fig, self.ax = plt.subplots(figsize=(6, 3.2), facecolor='#f5f7fa')
        self._sax()
        self.ax.text(0.5, 0.5, self.L['hint'], transform=self.ax.transAxes,
                     ha='center', va='center', color='#57606a', fontsize=10)
        self.cv = FigureCanvasTkAgg(self.fig, master=f)
        self.cv.get_tk_widget().pack(fill='both', expand=True, padx=4, pady=4)

    def _sax(self):
        self.ax.set_facecolor('#ffffff')
        self.ax.tick_params(colors='#57606a')
        for sp in self.ax.spines.values():
            sp.set_color('#d0d7de')
        self.ax.set_xlabel(self.L['time'], color='#57606a')
        self.ax.set_ylabel(self.L['level'], color='#57606a')

    def _draw_graph(self, t, sol, title=''):
        self.ax.clear(); self._sax()
        cols = ['#cf222e', '#1a7f37', '#0969da', '#bc4c00']
        labels = [self.L['dis_prog'], self.L['prot_lvl'],
                  self.L['immune'], self.L['drug_conc']]
        for i in range(sol.shape[1]):
            lw = 1.5 if i == 2 else 2
            ls = '--' if i == 2 else '-'
            self.ax.plot(t, sol[:, i], color=cols[i], lw=lw, ls=ls, label=labels[i])
        self.ax.set_ylim(0, 1.05)
        self.ax.grid(True, color='#eef1f5', lw=0.5)
        self.ax.legend(facecolor='#eef1f5', edgecolor='#d0d7de',
                       labelcolor='#1c2128', fontsize=8)
        if title:
            self.ax.set_title(title, color='#57606a', fontsize=9)
        self.cv.draw()

    def _build_3d(self, p):
        f = self._frm(p, self.L['visual_3d']); f.pack(fill='both', expand=True)
        # White medical background; narrow canvas centered (image is tall & narrow)
        holder = tk.Frame(f, bg=Theme.c('panel'))
        holder.pack(fill='both', expand=True, padx=4, pady=4)
        self.c3 = tk.Canvas(holder, bg=Theme.c('panel'), width=300, height=440,
                            highlightthickness=0)
        self.c3.pack()  # centered, not stretched
        self.human = AnatomyImageBody(self.c3)
        self.after(150, lambda: self.human.draw_default(self.L['hint']))



    def _build_controls(self, p):
        f = self._frm(p, self.L['mode']); f.pack(side='top', fill='x', pady=(0, 4))
        nb = ttk.Notebook(f); nb.pack(fill='x', padx=4, pady=4)
        st = ttk.Style()
        st.configure('TNotebook', background='#ffffff', borderwidth=0)
        st.configure('TNotebook.Tab', background='#eef1f5',
                      foreground=Theme.c('text_dim'), padding=[10, 4])
        st.map('TNotebook.Tab',
               background=[('selected', '#1a7f37')],
               foreground=[('selected', 'white')])

        is_ru = (self.lang == 'ru')
        mol = tk.Frame(nb, bg=Theme.c('panel')); nb.add(mol, text=self.L['mol_bio'])
        self._mol_tab(mol)
        # MolBio action buttons live INSIDE the MolBio tab, so they appear only
        # when MolBio is selected and disappear when switching to SynBio.
        mol_btns = tk.Frame(mol, bg=Theme.c('panel'))
        mol_btns.pack(fill='x', padx=4, pady=(4, 2))
        row_ab = tk.Frame(mol_btns, bg=Theme.c('panel')); row_ab.pack(fill='x')
        tk.Button(row_ab, text=self.L['analyze'], bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                  relief='solid', bd=1, font=('Arial', 9, 'bold'), cursor='hand2',
                  pady=3, command=self._run_api).pack(side='left', fill='x', expand=True, padx=(0, 2))
        tk.Button(row_ab, text=self.L['screen'], bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                  relief='solid', bd=1, font=('Arial', 9, 'bold'), cursor='hand2',
                  pady=3, command=self._run_screen).pack(side='left', fill='x', expand=True, padx=(2, 0))
        tk.Button(mol_btns, text='MolBio Workspace', bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                   relief='solid', bd=1, font=('Arial', 9, 'bold'), cursor='hand2',
                   pady=3, command=self._open_molbio
                   ).pack(fill='x', pady=(2, 0))

        syn = tk.Frame(nb, bg=Theme.c('panel')); nb.add(syn, text=self.L['syn_bio'])
        self._syn_tab(syn)

        dr = tk.Frame(f, bg=Theme.c('panel')); dr.pack(fill='x', padx=6, pady=(0, 4))
        tk.Label(dr, text=self.L['days'] + ':', bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                 font=('Arial', 9)).pack(side='left')
        self.dlbl = tk.Label(dr, text='100', bg=Theme.c('panel'), fg='#0969da',
                              font=('Arial', 10, 'bold'), width=4)
        self.dlbl.pack(side='right')
        tk.Scale(dr, from_=30, to=365, orient='horizontal', variable=self.days_var,
                 bg=Theme.c('panel'), fg=Theme.c('text'), highlightthickness=0,
                 troughcolor='#eef1f5', activebackground='#1a7f37',
                 command=lambda v: self.dlbl.config(text=v)
                 ).pack(side='left', fill='x', expand=True)

        bf = tk.Frame(f, bg=Theme.c('panel')); bf.pack(fill='x', padx=4, pady=(0, 6))
        tk.Button(bf, text=self.L['run'], bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                  relief='flat', font=('Arial', 12, 'bold'), cursor='hand2',
                  pady=7, command=self._run).pack(fill='x', pady=(0, 4))
        # Row: Reset + All Tools (moved here from top bar)
        row2 = tk.Frame(bf, bg=Theme.c('panel')); row2.pack(fill='x')
        tk.Button(row2, text=self.L['reset'], bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                  relief='solid', bd=1, font=('Arial', 10), cursor='hand2',
                  command=self._reset).pack(side='left', fill='x', expand=True, padx=(0, 2))
        tk.Button(row2, text=('Все инструменты' if self.lang == 'ru' else 'All tools'),
                  bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                  relief='solid', bd=1, font=('Arial', 10, 'bold'), cursor='hand2',
                  command=self._open_all_tools).pack(side='left', fill='x', expand=True, padx=(2, 0))

    def _mol_tab(self, p):
        tk.Label(p, text=self.L['drug_lbl'], bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                 font=('Arial', 10)).pack(anchor='w', padx=6, pady=(6, 2))
        self.drug_menu = ttk.Combobox(p, textvariable=self.drug_var,
                                       state='readonly', font=('Arial', 10))
        self.drug_menu.pack(fill='x', padx=6, pady=(0, 4))
        self.drug_menu['values'] = [self.L['sel_dis']]
        self.drug_menu.current(0)

        # Checkbox: disable preset drugs (so created drugs don't conflict)
        self.disable_preset_var = tk.IntVar(value=0)
        tk.Checkbutton(p, text=('Не использовать готовые лекарства'
                                  if self.lang == 'ru' else 'Disable preset drugs'),
                        variable=self.disable_preset_var,
                        bg=Theme.c('panel'), fg='#9a6700', selectcolor='#eef1f5',
                        activebackground='#ffffff', activeforeground='#9a6700',
                        font=('Arial', 8), command=self._toggle_preset
                        ).pack(anchor='w', padx=6, pady=(0, 4))

        tk.Label(p, text=self.L['dose_pct'], bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                 font=('Arial', 10)).pack(anchor='w', padx=6, pady=(0, 2))
        dr = tk.Frame(p, bg=Theme.c('panel')); dr.pack(fill='x', padx=6, pady=(0, 4))
        self.dose_lbl = tk.Label(dr, text='50', bg=Theme.c('panel'), fg='#0969da',
                                  font=('Arial', 11, 'bold'), width=3)
        self.dose_lbl.pack(side='right')
        tk.Scale(dr, from_=0, to=100, orient='horizontal', variable=self.drug_dose,
                 bg=Theme.c('panel'), fg=Theme.c('text'), highlightthickness=0,
                 troughcolor='#eef1f5', activebackground='#1a7f37',
                 command=lambda v: self.dose_lbl.config(text=v)
                 ).pack(side='left', fill='x', expand=True)

        # SMILES input + dose mg (compact row)
        tk.Label(p, text='SMILES:', bg=Theme.c('panel'), fg='#bc4c00',
                 font=('Arial', 9, 'bold')).pack(anchor='w', padx=6, pady=(4, 0))
        self.custom_drug = tk.StringVar()
        tk.Entry(p, textvariable=self.custom_drug, bg=Theme.c('panel2'), fg=Theme.c('text'),
                 insertbackground='#1c2128', relief='flat', font=('Arial', 9),
                 highlightthickness=1, highlightbackground=Theme.c('border')
                 ).pack(fill='x', padx=6, pady=(0, 4), ipady=3)

        dmr = tk.Frame(p, bg=Theme.c('panel')); dmr.pack(fill='x', padx=6, pady=(0, 4))
        tk.Label(dmr, text=self.L['dose_mg'] + ':', bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                 font=('Arial', 9)).pack(side='left')
        self.dmg_lbl = tk.Label(dmr, text='100', bg=Theme.c('panel'), fg='#0969da',
                                  font=('Arial', 10, 'bold'), width=5)
        self.dmg_lbl.pack(side='right')
        tk.Scale(dmr, from_=1, to=2000, orient='horizontal',
                 variable=self.dose_mg_var, bg=Theme.c('panel'), fg=Theme.c('text'),
                 highlightthickness=0, troughcolor='#eef1f5',
                 activebackground='#1a7f37',
                 command=lambda v: self.dmg_lbl.config(text=v)
                 ).pack(side='left', fill='x', expand=True)

        # Weight/height row
        wh = tk.Frame(p, bg=Theme.c('panel')); wh.pack(fill='x', padx=6, pady=(0, 2))
        tk.Label(wh, text=self.L['weight'], bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                 font=('Arial', 8)).pack(side='left')
        tk.Entry(wh, textvariable=self.weight_var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                 width=4, relief='flat', font=('Arial', 9),
                 highlightthickness=1, highlightbackground=Theme.c('border')
                 ).pack(side='left', padx=4)
        tk.Label(wh, text=self.L['height'], bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                 font=('Arial', 8)).pack(side='left', padx=(8, 0))
        tk.Entry(wh, textvariable=self.height_var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                 width=4, relief='flat', font=('Arial', 9),
                 highlightthickness=1, highlightbackground=Theme.c('border')
                 ).pack(side='left', padx=4)

        # Checkbox: analyze WITHOUT internet (RDKit only, for created molecules)
        self.force_rdkit_var = tk.IntVar(value=0)
        tk.Checkbutton(p, text=('Без интернета (RDKit) — для своих молекул'
                                  if self.lang == 'ru' else
                                  'Offline (RDKit only) — for created molecules'),
                        variable=self.force_rdkit_var,
                        bg=Theme.c('panel'), fg='#1a7f37', selectcolor='#eef1f5',
                        activebackground='#ffffff', activeforeground='#1a7f37',
                        font=('Arial', 8)
                        ).pack(anchor='w', padx=6, pady=(0, 4))

    def _syn_tab(self, p):
        # Synthetic Biology launcher tab — opens dedicated SynBio window
        # The old "mutation" controls (mut_menu, mut_val radios) are still needed
        # because _run() and _run_screen() reference them. We keep them hidden
        # but functional.
        self.mut_menu = ttk.Combobox(p, textvariable=self.mut_gene,
                                      state='readonly', font=('Arial', 10))
        # NOT packed (hidden) — but values still set for _run() compatibility
        self.mut_menu['values'] = HEALTHY_GENES
        self.mut_menu.current(0)
        self.mut_gene.set(HEALTHY_GENES[0])

        # Header
        is_ru = (self.lang == 'ru')
        title_text = 'Синтетическая биология' if is_ru else 'Synthetic Biology'
        desc_text = ('Создание генных цепей, CRISPR редактирование,\n'
                     'плазмиды, метаболическая инженерия.\n'
                     'Откроется отдельное окно — главное окно\n'
                     'продолжит работать.') if is_ru else (
                     'Gene circuits, CRISPR editing,\n'
                     'plasmids, metabolic engineering.\n'
                     'Opens in a separate window — main\n'
                     'window keeps working.')

        tk.Label(p, text=title_text, font=('Arial', 13, 'bold'),
                  bg=Theme.c('panel'), fg='#8250df').pack(anchor='w', padx=10, pady=(12, 4))
        tk.Label(p, text=desc_text, font=('Arial', 9),
                  bg=Theme.c('panel'), fg=Theme.c('text_dim'), justify='left'
                  ).pack(anchor='w', padx=10, pady=(0, 12))

        # Feature list
        features = (['Генные цепи (Repressilator, Toggle Switch и др.)',
                      'CRISPR/Cas9 редактирование генов',
                      'Дизайн плазмид (промоторы, гены, терминаторы)',
                      'Метаболическая инженерия']
                     if is_ru else
                     ['Gene circuits (Repressilator, Toggle Switch, etc.)',
                      'CRISPR/Cas9 gene editing',
                      'Plasmid design (promoters, genes, terminators)',
                      'Metabolic engineering'])
        for feat in features:
            tk.Label(p, text='  • ' + feat, font=('Arial', 9),
                      bg=Theme.c('panel'), fg=Theme.c('text'), anchor='w'
                      ).pack(anchor='w', padx=14, pady=1)

        # Big launch button
        start_text = 'НАЧАТЬ РАБОТУ' if is_ru else 'START'
        tk.Button(p, text=start_text, bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                   relief='solid', bd=1, font=('Arial', 13, 'bold'), cursor='hand2',
                   pady=10, command=self._open_synbio
                   ).pack(fill='x', padx=10, pady=(16, 8))

        # Mutation controls (kept for backward compatibility, in a small section)
        sep = tk.Frame(p, bg='#d0d7de', height=1)
        sep.pack(fill='x', padx=10, pady=8)

        legacy_title = 'Симуляция мутации (быстро)' if is_ru else 'Quick mutation simulation'
        tk.Label(p, text=legacy_title, font=('Arial', 10, 'bold'),
                  bg=Theme.c('panel'), fg='#1a7f37').pack(anchor='w', padx=10, pady=(4, 2))
        gene_lbl = 'Ген:' if is_ru else 'Gene:'
        tk.Label(p, text=gene_lbl, bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                  font=('Arial', 9)).pack(anchor='w', padx=10, pady=(2, 0))
        # Re-pack the combobox now (visible)
        self.mut_menu.pack(fill='x', padx=10, pady=(0, 4))

        mut_lbl = 'Мутация:' if is_ru else 'Mutation:'
        tk.Label(p, text=mut_lbl, bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                  font=('Arial', 9)).pack(anchor='w', padx=10, pady=(0, 2))
        row = tk.Frame(p, bg=Theme.c('panel')); row.pack(fill='x', padx=10, pady=(0, 8))
        tk.Radiobutton(row, text='OFF', variable=self.mut_val, value=0,
                       bg=Theme.c('panel'), fg='#1a7f37', selectcolor='#eef1f5',
                       activebackground='#ffffff',
                       font=('Arial', 10)).pack(side='left')
        tk.Radiobutton(row, text='ON', variable=self.mut_val, value=1,
                       bg=Theme.c('panel'), fg='#bc4c00', selectcolor='#eef1f5',
                       activebackground='#ffffff',
                       font=('Arial', 10, 'bold')).pack(side='left', padx=10)

    def _toggle_panel(self, panel_key):
        """Show or hide a main panel from the View menu (PanedWindow version)."""
        panel = self._panels.get(panel_key)
        if not panel:
            return
        pw = self._main_container
        if self._panel_vars[panel_key].get():
            # Show — re-add in the correct order position
            try:
                # Find where to insert based on panel_order
                visible = [k for k in self._panel_order
                           if k != panel_key and str(self._panels[k]) in
                           [str(p) for p in pw.panes()]]
                idx = self._panel_order.index(panel_key)
                # Insert after the last visible panel that comes before it
                before_panel = None
                for k in self._panel_order[idx+1:]:
                    if k in visible:
                        before_panel = self._panels[k]
                        break
                opts = {'minsize': 200}
                if panel_key == 'center':
                    opts['stretch'] = 'always'
                else:
                    opts['stretch'] = 'never'
                    opts['width'] = 290 if panel_key == 'left' else 360
                if before_panel is not None:
                    pw.add(panel, before=before_panel, **opts)
                else:
                    pw.add(panel, **opts)
            except Exception as _exc:
                biosim_log('warn', '_toggle_panel: {}'.format(_exc))
                pw.add(panel)
        else:
            try:
                pw.forget(panel)
            except Exception as _exc:
                biosim_log('warn', '_toggle_panel: {}'.format(_exc))

    def _toggle_preset(self):
        """Toggle preset drug list — when ON, clear and disable it."""
        if self.disable_preset_var.get() == 1:
            # Disable preset drugs
            self.drug_menu['values'] = ['-- disabled --']
            self.drug_menu.current(0)
            self.drug_var.set('')
            self.drug_menu.config(state='disabled')
        else:
            # Re-enable; restore based on current disease
            self.drug_menu.config(state='readonly')
            if self.cur_dis and self.cur_dis in DB:
                self._update_menus(self.cur_dis)
            else:
                self.drug_menu['values'] = [self.L['sel_dis']]
                self.drug_menu.current(0)

    def _update_menus(self, key):
        d = DB.get(key, {})
        drugs = list(d.get('drugs', {}).keys())
        # Respect disable_preset toggle
        if hasattr(self, 'disable_preset_var') and self.disable_preset_var.get() == 1:
            self.drug_menu['values'] = ['-- disabled --']
            self.drug_menu.current(0)
            self.drug_var.set('')
        else:
            self.drug_menu['values'] = drugs if drugs else ['--']
            if drugs:
                self.drug_menu.current(0); self.drug_var.set(drugs[0])
        dg = d.get('genes', []) + [g for g in HEALTHY_GENES if g not in d.get('genes', [])]
        self.mut_menu['values'] = dg
        if dg:
            self.mut_menu.current(0); self.mut_gene.set(dg[0])

    def _build_text(self, p):
        f = self._frm(p, self.L['text_res'])
        f.pack(side='bottom', fill='x', pady=(4, 0))
        self.rt = scrolledtext.ScrolledText(
            f, bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 10), relief='flat',
            wrap='word', highlightthickness=0, height=7)
        self.rt.pack(fill='x', padx=4, pady=4)
        for tag, fg, fnt in [
            ('ttl', '#0969da', ('Arial', 11, 'bold')),
            ('ok', '#1a7f37', ('Arial', 10)),
            ('wrn', '#9a6700', ('Arial', 10)),
            ('bad', '#cf222e', ('Arial', 10)),
            ('se', '#bc4c00', ('Arial', 10)),
            ('nrm', '#1c2128', ('Arial', 10)),
            ('dim', '#57606a', ('Arial', 9)),
            ('pub', '#0550ae', ('Arial', 9)),
            ('gene', '#1a7f37', ('Arial', 10)),
            ('lbl', '#57606a', ('Arial', 10, 'bold')),
        ]:
            self.rt.tag_config(tag, foreground=fg, font=fnt)
        # Make read-only but still selectable & copyable:
        # block any key that would edit, but allow copy/select shortcuts.
        self._make_readonly_copyable(self.rt)

    def _make_readonly_copyable(self, widget):
        """Allow selecting and copying text, but block all editing."""
        def on_key(event):
            # Allow copy (Ctrl/Cmd+C), select-all (Ctrl/Cmd+A), and navigation
            ctrl = (event.state & 0x4) != 0      # Control
            cmd = (event.state & 0x8) != 0 or (event.state & 0x10) != 0  # Command (Mac)
            if (ctrl or cmd) and event.keysym.lower() in ('c', 'a', 'insert'):
                return  # allow
            # Allow arrow keys, Home/End, PageUp/Down for navigation
            if event.keysym in ('Left', 'Right', 'Up', 'Down', 'Home', 'End',
                                  'Prior', 'Next', 'Shift_L', 'Shift_R',
                                  'Control_L', 'Control_R', 'Meta_L', 'Meta_R'):
                return  # allow
            return 'break'  # block everything else (typing, delete, backspace, paste)
        widget.bind('<Key>', on_key)
        # Block paste and cut explicitly
        widget.bind('<<Paste>>', lambda e: 'break')
        widget.bind('<<Cut>>', lambda e: 'break')
        # Block middle-click paste (X11)
        widget.bind('<Button-2>', lambda e: 'break')

    def _rw(self, lines):
        self.rt.delete('1.0', 'end')
        for txt, tag in lines:
            self.rt.insert('end', txt, tag)

    def _run(self):
        gene = self.mut_gene.get()
        drug = self.drug_var.get()
        try:
            mv = int(float(self.mut_val.get()))
            dose = int(float(self.drug_dose.get()))
            days = int(float(self.days_var.get()))
        except (ValueError, TypeError) as _exc:
            biosim_log('warn', '_run: {}'.format(_exc))
            mv, dose, days = 0, 50, 100
        d = self.cur_dat

        if d and 'params' in d:
            t, sol = sim_disease(d, drug, dose, gene, mv, days)
            fD = float(sol[-1][0]); fP = float(sol[-1][1]); fI = float(sol[-1][2])
            name = d.get('ru', '') if self.lang == 'ru' else d.get('en', '')
            self._draw_graph(t, sol, name)
            self.human.draw_disease(d, fD, self.lang)
            self._write_preset(d, name, fD, fP, fI, drug, dose, gene, mv)
        elif mv == 1 and gene:
            t, sol = sim_mutation(gene, days)
            fD = float(sol[-1][0]); fP = float(sol[-1][1]); fI = float(sol[-1][2])
            self._draw_graph(t, sol, self.L['mut_only'] + ': ' + gene)
            self.human.draw_mutation(gene, fD)
            self._write_mutation(gene, fD, fP, fI)
        else:
            messagebox.showinfo('VIREON', self.L['hint'])

    def _write_preset(self, d, name, D, P, I, drug, dose, gene, mv):
        st = 'ok' if D < 0.3 else ('wrn' if D < 0.6 else 'bad')
        sl = self.L['low'] if D < 0.3 else (self.L['moderate'] if D < 0.6 else self.L['high'])
        lines = [
            (self.L['results'] + '\n', 'ttl'),
            ('-' * 26 + '\n', 'dim'),
            (self.L['dis_prog'] + ': ', 'nrm'),
            ('{:.2f}  ({})\n'.format(D, sl), st),
            (self.L['prot_lvl'] + ': {:.2f}\n'.format(P), 'nrm'),
            (self.L['immune'] + ': {:.2f}\n\n'.format(I), 'nrm'),
        ]
        if drug and drug in d.get('drugs', {}):
            dd = d['drugs'][drug]
            eff = dd['eff'] * dose / 100.0
            et = 'ok' if eff > 0.4 else 'wrn'
            lines += [(drug + ' ({}%)\n'.format(dose), 'ttl'),
                      ('  ' + self.L['drug_eff'] + ': {:.0%}\n'.format(eff), et)]
            se_key = 'se_ru' if self.lang == 'ru' else 'se_en'
            ses = dd.get(se_key, [])
            if ses:
                lines.append(('\n' + self.L['side_eff'] + ':\n', 'se'))
                for s in ses:
                    lines.append(('  - ' + s + '\n', 'se'))
        if mv == 1 and gene:
            ek = 'ru' if self.lang == 'ru' else 'en'
            et = GENE_INFO.get(gene, {}).get(ek, '')
            lines += [('\n' + self.L['mutation'] + ': ' + gene + ' ON\n', 'bad'),
                      ('  ' + et + '\n', 'bad'),
                      ('  ' + self.L['prog25'] + '\n', 'bad')]
        key = self.cur_dis
        if key and key in DB:
            ld = DB[key]
            lines += [('\n' + self.L['survival'] + ': {}%\n'.format(ld['surv']), 'ok'),
                      (self.L['recovery'] + ': {}%\n'.format(ld['recov']), 'ok'),
                      (self.L['management'] + ': {}%\n'.format(ld['mgmt']), 'ok')]
        self._rw(lines)

    def _write_mutation(self, gene, D, P, I):
        ek = 'ru' if self.lang == 'ru' else 'en'
        gi = GENE_INFO.get(gene, {})
        et = gi.get(ek, '')
        risk = float(gi.get('risk', 0.5))
        st = 'ok' if D < 0.25 else ('wrn' if D < 0.5 else 'bad')
        lines = [
            (self.L['mut_only'] + '\n', 'ttl'),
            ('-' * 26 + '\n', 'dim'),
            (self.L['genes'] + ': ' + gene + '\n', 'bad'),
            ('  ' + et + '\n\n', 'nrm'),
            (self.L['risk'] + ': {:.0%}\n'.format(risk), st),
            (self.L['prot_lvl'] + ': {:.2f}\n'.format(P), 'nrm'),
            (self.L['immune'] + ': {:.2f}\n'.format(I), 'nrm'),
        ]
        self._rw(lines)

    def _run_screen(self):
        """Create new drug molecules for the current disease."""
        if not self.cur_dat or 'genes' not in self.cur_dat:
            messagebox.showinfo('VIREON', self.L['sel_dis'])
            return
        if not _ensure_rdkit():
            messagebox.showinfo('VIREON', 'RDKit не установлен. Установите: pip3.8 install rdkit')
            return
        genes = self.cur_dat['genes']
        self._rw([(self.L['screen'] + '...\n', 'wrn')])
        self.update()
        def worker():
            def cb(s, t, m):
                self.after(0, lambda: self._rw(
                    [(self.L['screen'] + ':\n', 'wrn'),
                     ('[{}/{}] {}\n'.format(s, t, m), 'dim')]))
            results = create_drugs_for_disease(genes, cb, disease_name=self.cur_dis or '')
            self.after(0, lambda: self._show_created_drugs(results))
        threading.Thread(target=worker, daemon=True).start()

    def _show_created_drugs(self, results):
        """Display newly created drug molecules."""
        if not results:
            self._rw([(self.L['not_found'] + '\n', 'dim')])
            return
        lines = [(self.L['screen'] + ' (de novo)\n', 'ttl'),
                  ('-' * 30 + '\n', 'dim'),
                  ('Generated via RDKit mutations + drug-likeness scoring\n', 'dim')]
        if _ensure_deepchem():
            lines.append(('DeepChem: ON  ', 'ok'))
        if _ensure_torch():
            lines.append(('PyTorch: ON  ', 'ok'))
        lines.append(('RDKit: ON\n\n', 'ok'))
        
        for idx, R in enumerate(results, 1):
            lines.append(('TARGET {}: {} (gene {})\n'.format(
                idx, R['target'][:50], R.get('gene', '?')), 'ttl'))
            lines.append(('  Function: {}\n'.format(R.get('protein_function', '')[:120]), 'dim'))
            if R.get('protein_acc'):
                lines.append(('  UniProt: {} | AlphaFold: alphafold.ebi.ac.uk/entry/{}\n'.format(
                    R['protein_acc'], R['protein_acc']), 'dim'))
            lines.append(('  Total generated: {} molecules, top 5:\n\n'.format(
                R.get('total_generated', 0)), 'dim'))
            
            for j, mol in enumerate(R['molecules'][:5], 1):
                score_tag = 'ok' if mol['score'] > 60 else ('wrn' if mol['score'] > 40 else 'bad')
                lines.append(('  CANDIDATE {}: Score {}/100\n'.format(j, mol['score']), score_tag))
                lines.append(('    SMILES: {}\n'.format(mol['smiles'][:70]), 'nrm'))
                lines.append(('    Formula: {}  MW: {} Da\n'.format(mol['formula'], mol['mw']), 'nrm'))
                lines.append(('    LogP: {}  TPSA: {}  QED: {}\n'.format(
                    mol['logp'], mol['tpsa'], mol['qed']), 'dim'))
                lines.append(('    HBD: {}  HBA: {}  RotBonds: {}  Lipinski violations: {}\n'.format(
                    mol['hbd'], mol['hba'], mol['rtb'], mol['lipinski_violations']), 'dim'))
                lines.append(('    Parent: {} ({})\n'.format(mol.get('parent_name', '?'), mol.get('mutation', '?')), 'dim'))
                for w in mol.get('tox_warnings', []):
                    lines.append(('    ! ' + w + '\n', 'bad'))
                lines.append(('\n', 'nrm'))
            lines.append(('\n', 'nrm'))
        
        lines.append(('To test a candidate: copy its SMILES and paste into "Drug analysis" field below.\n', 'dim'))
        self._rw(lines)


    def _run_api(self):
        query = self.custom_drug.get().strip()
        if not query:
            messagebox.showinfo('VIREON', self.L['drug_input'])
            return
        try:
            dose_mg = int(self.dose_mg_var.get())
            weight = int(self.weight_var.get())
            height = int(self.height_var.get())
        except (tk.TclError, ValueError) as _exc:
            biosim_log('warn', '_run_api: {}'.format(_exc))
            dose_mg = 100; weight = 70; height = 170
        force_rdkit = bool(self.force_rdkit_var.get())
        self._rw([(self.L['loading'] + '...\n', 'wrn')])
        self.update()
        def worker():
            def cb(s, t, m):
                self.after(0, lambda: self._rw(
                    [(self.L['loading'] + ':\n', 'wrn'),
                     ('[{}/{}] {}\n'.format(s, t, m), 'dim')]))
            R = analyze_drug(query, self.cur_dat, dose_mg, weight, height, cb, force_rdkit=force_rdkit)
            self.after(0, lambda: self._show_api(R))
        threading.Thread(target=worker, daemon=True).start()

    def _show_api(self, R):
        if R['error']:
            self._rw([(self.L['not_found'] + '\n', 'bad'), (R['error'] + '\n', 'dim')])
            return
        self._api_result = R

        if self.cur_dat and 'params' in self.cur_dat and R['ode_params']:
            d = self.cur_dat
            try:
                dose_frac = int(float(self.drug_dose.get())) / 100.0
                days = int(float(self.days_var.get()))
            except (ValueError, TypeError) as _exc:
                biosim_log('warn', '_show_api: {}'.format(_exc))
                dose_frac, days = 0.5, 100
            t, sol = sim_disease_api(d, R['ode_params'], dose_frac, days)
            fD = float(sol[-1][0])
            self._draw_graph(t, sol, R['name'])
            self.human.draw_disease(d, fD, self.lang)

        if R['alphafold']:
            self.after(100, lambda: self.human.draw_alphafold(R['alphafold']))

        lines = [
            (R['name'] + '\n', 'ttl'),
            ('-' * 26 + '\n', 'dim'),
        ]
        if R.get('source') == 'rdkit':
            lines.append(('Source: RDKit (novel molecule)\n', 'wrn'))
        lines.append(('Formula: ' + R['formula'] + '  MW: {:.1f} Da\n'.format(R['mw']), 'nrm'))
        lines.append(('SMILES: ' + R['smiles'][:55] + '\n', 'dim'))
        if R['desc']:
            lines.append((R['desc'][:280] + '\n\n', 'nrm'))

        li = R['lipinski']
        lines.append((self.L['lipinski'] + ':\n', 'ttl'))
        for rule_name, passed, val in li.get('rules', []):
            lines.append(('  {} {}\n'.format('OK' if passed else 'FAIL', rule_name),
                          'ok' if passed else 'bad'))
        lines.append(('  QED: {:.3f}  Oral: {}\n\n'.format(
            li.get('score', 0), 'OK' if li.get('oral_ok') else 'Poor'),
            'ok' if li.get('oral_ok') else 'wrn'))

        bio = R.get('biophysics', {})
        if bio:
            lines.append((self.L['biophysics'] + ':\n', 'ttl'))
            if bio.get('bioavailability') is not None:
                lines.append(('  ' + self.L['bioavail'] + ': {:.1f}%\n'.format(bio['bioavailability']),
                              'ok' if bio['bioavailability'] > 50 else 'wrn'))
            if bio.get('half_life_h') is not None:
                lines.append(('  ' + self.L['half_life'] + ': {} h\n'.format(bio['half_life_h']), 'nrm'))
            if bio.get('vd_total_L') is not None:
                lines.append(('  ' + self.L['vd'] + ': {} L\n'.format(bio['vd_total_L']), 'nrm'))
            if bio.get('binding_dG') is not None:
                lines.append(('  ' + self.L['binding'] + ': {} kcal/mol\n'.format(bio['binding_dG']), 'ok'))
            if bio.get('affinity_nM'):
                lines.append(('  Affinity (Kd): {} nM\n'.format(bio['affinity_nM']), 'ok'))
            lines.append(('\n', 'nrm'))

        if R['targets']:
            lines.append((self.L['targets'] + ':\n', 'ttl'))
            for tgt in R['targets'][:5]:
                ic = tgt.get('ic50_nm')
                ics = '  IC50={} nM'.format(ic) if ic else ''
                lines.append(('  {} (pChEMBL={:.1f}{})\n'.format(
                    tgt['name'], tgt['pchembl'], ics), 'nrm'))
            lines.append(('\n', 'nrm'))

        if R['alphafold']:
            lines.append((self.L['alphafold'] + ':\n', 'ttl'))
            for af in R['alphafold']:
                conf = float(af.get('confidence', 0))
                tag = 'ok' if conf > 70 else ('wrn' if conf > 50 else 'bad')
                lines.append(('  {} {} conf={:.0f}%\n'.format(
                    af.get('gene', ''), af.get('protein', ''), conf), tag))
            lines.append(('\n', 'nrm'))
            # Show 3D buttons for each protein
            self._pending_3d_buttons = R['alphafold']

        if R['gene_data']:
            lines.append((self.L['dna'] + ':\n', 'ttl'))
            for gene, gd in R['gene_data'].items():
                lines.append(('  ' + gene + ': ' + gd.get('fullname', '')[:55] + '\n', 'nrm'))
                if gd.get('chrom'):
                    lines.append(('  Chr:{} Loc:{} GeneID:{}\n'.format(
                        gd['chrom'], gd.get('loc', ''), gd.get('id', '')), 'dim'))
            lines.append(('\n', 'nrm'))
            self._pending_dna_buttons = R['gene_data']

        admet = R['admet']
        if admet:
            lines.append(('ADMET:\n', 'ttl'))
            if admet.get('max_phase') is not None:
                lines.append(('  ' + self.L['phase'] + ': {}\n'.format(admet['max_phase']),
                              'ok' if admet['max_phase'] == 4 else 'nrm'))
            if admet.get('oral') is not None:
                lines.append(('  ' + self.L['oral'] + ': {}\n'.format(
                    'Yes' if admet['oral'] else 'No'),
                    'ok' if admet['oral'] else 'wrn'))
            lines.append(('\n', 'nrm'))

        ov = R['overdose']
        if ov:
            lt = 'ok' if ov['level'] == 'low' else ('wrn' if ov['level'] == 'moderate' else 'bad')
            lines += [
                (self.L['overdose'] + ':\n', 'ttl'),
                ('  Weight: {} kg, Height: {} cm, BMI: {}\n'.format(
                    ov['weight'], ov['height'], ov['bmi']), 'dim'),
                ('  BSA: {} m^2\n'.format(ov['bsa']), 'dim'),
                ('  ' + self.L['ti'] + ': {}\n'.format(ov['ti']), lt),
                ('  Risk: {}\n'.format(ov['level']), lt),
                ('  Lethal estimate: {} mg total ({} mg/kg)\n'.format(
                    ov['lethal_total'], ov['lethal_per_kg']), 'nrm'),
            ]
            if ov['bb']:
                lines.append(('  ' + self.L['blackbox'] + '\n', 'bad'))
            for sym in ov.get('symptoms', []):
                lines.append(('  - ' + sym + '\n', 'se'))
            lines.append(('\n', 'nrm'))

        out = R['outcomes']
        if out:
            lines += [
                ('Outcomes:\n', 'ttl'),
                ('  ' + self.L['survival'] + ': {}%\n'.format(out['surv']),
                 'ok' if out['surv'] > 60 else 'wrn'),
                ('  ' + self.L['recovery'] + ': {}%\n'.format(out['recov']),
                 'ok' if out['recov'] > 30 else 'wrn'),
                ('  ' + self.L['management'] + ': {}%\n'.format(out['mgmt']),
                 'ok' if out['mgmt'] > 50 else 'wrn'),
                ('\n', 'nrm'),
            ]

        if R['pathways']:
            lines.append((self.L['pathways'] + ':\n', 'ttl'))
            for pw in R['pathways'][:6]:
                lines.append(('  - ' + pw[:110] + '\n', 'nrm'))
            lines.append(('\n', 'nrm'))

        if R['publications']:
            lines.append((self.L['pubs'] + ':\n', 'ttl'))
            for pub in R['publications']:
                lines.append(('  [{}] {}\n'.format(pub['year'], pub['title'][:95]), 'pub'))
            lines.append(('\n', 'nrm'))

        self._rw(lines)
        # Show 3D action buttons under the text
        self._show_3d_buttons(R)

    def _show_3d_buttons(self, R):
        """Show buttons to open 2.5D viewers for proteins, DNA, molecule."""
        # Clear any previous buttons frame
        if hasattr(self, '_btns_frame') and self._btns_frame:
            try:
                self._btns_frame.destroy()
            except Exception as _exc:
                biosim_log('warn', '_show_3d_buttons: {}'.format(_exc))
            self._btns_frame = None
        # Find parent of self.rt
        parent = self.rt.master
        self._btns_frame = tk.Frame(parent, bg=Theme.c('panel'))
        self._btns_frame.pack(fill='x', padx=4, pady=(0, 4))
        tk.Label(self._btns_frame, text='2.5D Viewers:', bg=Theme.c('panel'),
                  fg='#0969da', font=('Arial', 9, 'bold')).pack(anchor='w', padx=4, pady=(2, 2))
        row = tk.Frame(self._btns_frame, bg=Theme.c('panel'))
        row.pack(fill='x', padx=4, pady=(0, 4))
        # Protein buttons
        for af in R.get('alphafold', [])[:3]:
            acc = af.get('acc', '')
            gene = af.get('gene', '')
            pname = af.get('protein', '')
            if acc:
                label = '3D: ' + (gene or pname[:10])
                tk.Button(row, text=label, bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                          font=('Arial', 8, 'bold'), cursor='hand2', padx=6,
                          command=lambda a=acc, n=pname, g=gene: Protein3DViewer(self, a, n, g)
                          ).pack(side='left', padx=2)
        # DNA buttons
        for gene_sym, gd in list(R.get('gene_data', {}).items())[:3]:
            tk.Button(row, text='DNA: ' + gene_sym, bg='#1a7f37', highlightbackground='#1a7f37', fg='#f5f7fa', relief='flat',
                      font=('Arial', 8, 'bold'), cursor='hand2', padx=6,
                      command=lambda g=gene_sym, d=gd: DNAViewer(self, g, d)
                      ).pack(side='left', padx=2)

    # ═══════════════════════════════════════════════════
    # MAIN WINDOW ENHANCEMENTS
    # Disease comparison, Treatment timeline, Lab values,
    # DICOM viewer, AI Doctor, Anatomy quiz
    # ═══════════════════════════════════════════════════

    def _open_disease_comparison(self):
        """Compare 2-3 diseases side-by-side."""
        # Get available diseases
        diseases = dict(DB)
        diseases.update(self._custom_diseases)

        # Dialog to select 2-3 diseases
        sel_dlg = tk.Toplevel(self)
        sel_dlg.title('' + ('Сравнение болезней' if self.lang == 'ru'
                                  else 'Disease Comparison'))
        sel_dlg.geometry('500x500')
        sel_dlg.configure(bg=Theme.c('bg'))
        sel_dlg.transient(self)

        tk.Label(sel_dlg, text=('Выбери 2-3 болезни для сравнения'
                                          if self.lang == 'ru' else
                                          'Choose 2-3 diseases to compare'),
                  bg=Theme.c('bg'), fg='#0969da',
                  font=('Arial', 12, 'bold')).pack(pady=10)

        # Scrollable checkbox list
        list_frame = tk.Frame(sel_dlg, bg=Theme.c('bg'))
        list_frame.pack(fill='both', expand=True, padx=20, pady=4)
        canv = tk.Canvas(list_frame, bg=Theme.c('bg'), highlightthickness=0)
        sb = tk.Scrollbar(list_frame, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True)
        sb.pack(side='right', fill='y')
        inner = tk.Frame(canv, bg=Theme.c('bg'))
        inner.bind('<Configure>',
                    lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=inner, anchor='nw', width=440)

        dis_vars = {}
        for key, dis in sorted(diseases.items()):
            v = tk.IntVar(value=0)
            dis_vars[key] = v
            name = dis.get(self.lang) or dis.get('ru') or key
            tk.Checkbutton(inner, text=name + '  (' + key + ')',
                            variable=v,
                            bg=Theme.c('bg'), fg=Theme.c('text'),
                            selectcolor='#eef1f5',
                            activebackground='#f5f7fa',
                            font=('Arial', 10), anchor='w'
                            ).pack(fill='x', padx=4, pady=1)

        def do_compare():
            selected = [k for k, v in dis_vars.items() if v.get()]
            if len(selected) < 2:
                return
            if len(selected) > 3:
                selected = selected[:3]
            sel_dlg.destroy()
            self._show_disease_comparison_window(selected, diseases)

        tk.Button(sel_dlg, text='▶ ' + ('Сравнить' if self.lang == 'ru'
                                            else 'Compare'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 11, 'bold'), cursor='hand2',
                   padx=16, pady=5, command=do_compare).pack(pady=10)

    def _show_disease_comparison_window(self, keys, diseases):
        """Show side-by-side comparison of selected diseases."""
        dlg = tk.Toplevel(self)
        dlg.title('Disease Comparison')
        dlg.geometry('{}x680'.format(360 * len(keys) + 60))
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Сравнение' if self.lang == 'ru'
                                      else 'Comparison') +
                        ': ' + ' · '.join(keys),
                  bg=Theme.c('bg'), fg='#0969da',
                  font=('Arial', 13, 'bold')).pack(pady=8)

        # Columns
        cols_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        cols_frame.pack(fill='both', expand=True, padx=10, pady=4)

        for key in keys:
            dis = diseases[key]
            col = tk.Frame(cols_frame, bg=Theme.c('panel'),
                              highlightthickness=2,
                              highlightbackground=Theme.c('border'))
            col.pack(side='left', fill='both', expand=True, padx=4, pady=4)

            # Header
            name = dis.get(self.lang) or dis.get('ru') or key
            tk.Label(col, text='' + name,
                      bg=Theme.c('panel'), fg='#0550ae',
                      font=('Arial', 12, 'bold'),
                      wraplength=320, justify='center'
                      ).pack(padx=10, pady=(10, 4))
            tk.Label(col, text=key,
                      bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                      font=('Courier', 9, 'italic')).pack()

            # Sections — use actual DB field names
            is_ru = (self.lang == 'ru')
            tk.Label(col, text=('Описание' if is_ru else 'Description'),
                      bg=Theme.c('panel'), fg='#1a7f37', font=('Arial', 10, 'bold')
                      ).pack(anchor='w', padx=12, pady=(8, 2))
            desc = dis.get('desc_ru' if is_ru else 'desc_en', dis.get('desc_ru', '?'))
            tk.Label(col, text=str(desc)[:250], bg=Theme.c('panel'), fg=Theme.c('text'),
                      font=('Arial', 9), wraplength=320, justify='left', anchor='w'
                      ).pack(anchor='w', padx=14, pady=(0, 2), fill='x')

            tk.Label(col, text=('Гены' if is_ru else 'Genes'),
                      bg=Theme.c('panel'), fg='#8250df', font=('Arial', 10, 'bold')
                      ).pack(anchor='w', padx=12, pady=(8, 2))
            genes = dis.get('genes', [])
            tk.Label(col, text=', '.join(genes) if genes else '—',
                      bg=Theme.c('panel'), fg=Theme.c('text'), font=('Arial', 9),
                      wraplength=320, justify='left', anchor='w'
                      ).pack(anchor='w', padx=14, pady=(0, 2), fill='x')

            tk.Label(col, text=('Поражённые органы' if is_ru else 'Affected organs'),
                      bg=Theme.c('panel'), fg='#bc4c00', font=('Arial', 10, 'bold')
                      ).pack(anchor='w', padx=12, pady=(8, 2))
            organs = dis.get('organs_ru' if is_ru else 'organs_en', dis.get('organs_ru', []))
            tk.Label(col, text=', '.join(organs) if organs else '—',
                      bg=Theme.c('panel'), fg=Theme.c('text'), font=('Arial', 9),
                      wraplength=320, justify='left', anchor='w'
                      ).pack(anchor='w', padx=14, pady=(0, 2), fill='x')

            tk.Label(col, text=('Лечение' if is_ru else 'Treatment'),
                      bg=Theme.c('panel'), fg='#9a6700', font=('Arial', 10, 'bold')
                      ).pack(anchor='w', padx=12, pady=(8, 2))
            drugs = dis.get('drugs', {})
            drug_names = ', '.join(drugs.keys()) if drugs else '—'
            tk.Label(col, text=drug_names, bg=Theme.c('panel'), fg=Theme.c('text'),
                      font=('Arial', 9), wraplength=320, justify='left', anchor='w'
                      ).pack(anchor='w', padx=14, pady=(0, 2), fill='x')

            # Survival/recovery stats
            stats = []
            if 'surv' in dis:
                stats.append(('Выживаемость' if is_ru else 'Survival') + ': {}%'.format(dis['surv']))
            if 'recov' in dis:
                stats.append(('Выздоровление' if is_ru else 'Recovery') + ': {}%'.format(dis['recov']))
            if stats:
                tk.Label(col, text=' · '.join(stats),
                          bg=Theme.c('panel'), fg='#cf222e',
                          font=('Courier', 9, 'bold'), wraplength=320
                          ).pack(anchor='w', padx=14, pady=6)

        # Close button
        tk.Button(dlg, text=('Закрыть' if self.lang == 'ru' else 'Close'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy
                   ).pack(fill='x', padx=14, pady=8)


    def _open_lab_values(self):
        """Show simulated lab values panel."""
        if not self.cur_dis:
            messagebox.showinfo('VIREON',
                ('Сначала выбери болезнь!' if self.lang == 'ru'
                  else 'Choose a disease first!'))
            return

        dlg = tk.Toplevel(self)
        dlg.title('' + ('Лабораторные показатели' if self.lang == 'ru'
                              else 'Lab Values'))
        dlg.geometry('800x620')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Лабораторные показатели' if self.lang == 'ru'
                                       else 'Lab Values'),
                  bg=Theme.c('bg'), fg='#1a7f37',
                  font=('Arial', 14, 'bold')).pack(pady=8)
        tk.Label(dlg, text=(self.cur_dat.get('name_' + self.lang) or
                                self.cur_dat.get('name_ru', '?')),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 10, 'italic')).pack()

        # Simulate lab values based on disease (using simple heuristics)
        import random
        rng = random.Random(hash(self.cur_dis or 'x'))

        # Reference values: (name, units, low_normal, high_normal)
        # Then simulate based on disease impact
        all_labs = [
            ('Glucose', 'mg/dL', 70, 100, 'metabolic'),
            ('HbA1c', '%', 4.0, 5.7, 'metabolic'),
            ('Insulin', 'μU/mL', 2.6, 24.9, 'metabolic'),
            ('Creatinine', 'mg/dL', 0.6, 1.2, 'kidney'),
            ('BUN', 'mg/dL', 7, 20, 'kidney'),
            ('eGFR', 'mL/min', 90, 120, 'kidney'),
            ('ALT', 'U/L', 7, 56, 'liver'),
            ('AST', 'U/L', 10, 40, 'liver'),
            ('Bilirubin', 'mg/dL', 0.1, 1.2, 'liver'),
            ('WBC', '×10³/μL', 4.5, 11, 'blood'),
            ('RBC', '×10⁶/μL', 4.5, 5.9, 'blood'),
            ('Hgb', 'g/dL', 13.5, 17.5, 'blood'),
            ('Platelets', '×10³/μL', 150, 400, 'blood'),
            ('Cholesterol', 'mg/dL', 0, 200, 'lipid'),
            ('LDL', 'mg/dL', 0, 100, 'lipid'),
            ('HDL', 'mg/dL', 40, 60, 'lipid'),
            ('Triglycerides', 'mg/dL', 0, 150, 'lipid'),
            ('CRP', 'mg/L', 0, 3, 'inflammation'),
            ('TSH', 'mIU/L', 0.4, 4.0, 'thyroid'),
        ]

        # Categorize organs of current disease
        affected_organs = self.cur_dat.get('organs', '').lower() if self.cur_dat else ''
        affected_categories = set()
        if 'pancrea' in affected_organs or 'diabet' in (self.cur_dis or '').lower():
            affected_categories.add('metabolic')
        if 'kidney' in affected_organs or 'почк' in affected_organs:
            affected_categories.add('kidney')
        if 'liver' in affected_organs or 'печен' in affected_organs:
            affected_categories.add('liver')
        if 'blood' in affected_organs or 'кров' in affected_organs:
            affected_categories.add('blood')
        if 'heart' in affected_organs or 'сердц' in affected_organs:
            affected_categories.add('lipid')

        # Display labs
        sb_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        sb_frame.pack(fill='both', expand=True, padx=14, pady=8)
        canv = tk.Canvas(sb_frame, bg=Theme.c('bg'), highlightthickness=0)
        sb = tk.Scrollbar(sb_frame, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True)
        sb.pack(side='right', fill='y')
        inner = tk.Frame(canv, bg=Theme.c('bg'))
        inner.bind('<Configure>',
                    lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=inner, anchor='nw', width=760)

        for name, units, low, high, category in all_labs:
            # Simulate value
            if category in affected_categories:
                # Abnormal — pick out of range
                if rng.random() < 0.5:
                    val = high * (1.2 + rng.random() * 0.8)  # high
                    status = 'HIGH'
                    color = '#cf222e'
                else:
                    val = low * (0.3 + rng.random() * 0.5)  # low
                    status = 'LOW'
                    color = '#0969da'
            else:
                # Normal range
                val = low + rng.random() * (high - low)
                status = 'OK'
                color = '#1a7f37'

            row = tk.Frame(inner, bg=Theme.c('panel'))
            row.pack(fill='x', padx=4, pady=1)

            tk.Label(row, text=name, bg=Theme.c('panel'), fg=Theme.c('text'),
                      font=('Arial', 10, 'bold'), width=14, anchor='w'
                      ).pack(side='left', padx=8, pady=4)
            tk.Label(row, text='{:.2f} {}'.format(val, units),
                      bg=Theme.c('panel'), fg=color,
                      font=('Courier', 11, 'bold'), width=14, anchor='w'
                      ).pack(side='left', padx=4)
            tk.Label(row,
                      text='ref: ' + '{:.1f}–{:.1f}'.format(low, high),
                      bg=Theme.c('panel'), fg=Theme.c('text_dim'),
                      font=('Courier', 9), width=18, anchor='w'
                      ).pack(side='left', padx=4)
            tk.Label(row, text=status, bg=Theme.c('panel'), fg=color,
                      font=('Arial', 9, 'bold'), width=8, anchor='w'
                      ).pack(side='left', padx=4)

        tk.Button(dlg, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy
                   ).pack(fill='x', padx=14, pady=8)

    def _open_dicom_viewer(self):
        """Simple medical image viewer (DICOM-like)."""
        try:
            from PIL import Image, ImageTk
        except ImportError:
            messagebox.showerror('VIREON', 'Pillow not installed')
            return

        from tkinter import filedialog
        filepath = filedialog.askopenfilename(
            title='Open medical image',
            filetypes=[('Images', '*.dcm *.png *.jpg *.jpeg *.tif *.tiff *.bmp'),
                        ('All files', '*.*')])
        if not filepath:
            return

        try:
            # Try opening as regular image first; DICOM requires pydicom
            if filepath.lower().endswith('.dcm'):
                try:
                    import pydicom
                    ds = pydicom.dcmread(filepath)
                    arr = ds.pixel_array
                    # Normalize
                    arr = ((arr - arr.min()) / max(1, arr.max() - arr.min()) * 255).astype('uint8')
                    img = Image.fromarray(arr)
                    metadata = '\n'.join([
                        'Patient: ' + str(ds.get('PatientName', '?')),
                        'Modality: ' + str(ds.get('Modality', '?')),
                        'Body part: ' + str(ds.get('BodyPartExamined', '?')),
                        'Study date: ' + str(ds.get('StudyDate', '?')),
                    ])
                except ImportError:
                    messagebox.showerror('VIREON',
                        'pydicom not installed. Install: pip install pydicom')
                    return
                except Exception as e:
                    messagebox.showerror('VIREON', 'DICOM error: ' + str(e))
                    return
            else:
                img = Image.open(filepath)
                metadata = 'Size: {}×{}\nMode: {}\nFormat: {}'.format(
                    img.size[0], img.size[1], img.mode, img.format or '?')
        except Exception as e:
            messagebox.showerror('VIREON', 'Cannot open: ' + str(e))
            return

        dlg = tk.Toplevel(self)
        dlg.title('Medical Image Viewer')
        dlg.geometry('960x720')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text='' + os.path.basename(filepath),
                  bg=Theme.c('bg'), fg='#0550ae',
                  font=('Arial', 13, 'bold')).pack(pady=8)

        main = tk.Frame(dlg, bg=Theme.c('bg'))
        main.pack(fill='both', expand=True, padx=10, pady=6)

        # Image
        img_frame = tk.Frame(main, bg=Theme.c('bg'))
        img_frame.pack(side='left', fill='both', expand=True, padx=4)

        # Resize to fit
        w, h = img.size
        max_dim = 600
        scale = min(max_dim / w, max_dim / h, 1.0)
        img_disp = img.resize((int(w * scale), int(h * scale)), Image.LANCZOS)
        photo = ImageTk.PhotoImage(img_disp)
        img_canvas = tk.Canvas(img_frame, bg=Theme.c('panel'),
                                   width=img_disp.size[0],
                                   height=img_disp.size[1],
                                   highlightthickness=2,
                                   highlightbackground=Theme.c('border'))
        img_canvas.pack(padx=4, pady=4)
        img_canvas.image = photo
        img_canvas.create_image(img_disp.size[0] // 2,
                                    img_disp.size[1] // 2,
                                    image=photo)

        # Metadata panel
        meta_frame = tk.LabelFrame(main,
                                       text='  ' + ('Метаданные' if self.lang == 'ru'
                                                        else 'Metadata') + '  ',
                                       bg=Theme.c('bg'), fg='#0550ae',
                                       font=('Arial', 10, 'bold'), relief='flat')
        meta_frame.pack(side='right', fill='y', padx=4)
        tk.Label(meta_frame, text=metadata,
                  bg=Theme.c('bg'), fg=Theme.c('text'),
                  font=('Courier', 9), anchor='w', justify='left'
                  ).pack(padx=10, pady=10, anchor='w')

        tk.Button(dlg, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy
                   ).pack(fill='x', padx=14, pady=8)







    def _open_patient_profile(self):
        """Patient profile dialog with medical history, allergies, medications."""
        # Load or initialize profile
        if not hasattr(self, '_patient_profile'):
            self._patient_profile = self._load_patient_profile()

        dlg = tk.Toplevel(self)
        dlg.title('' + ('Профиль пациента' if self.lang == 'ru' else 'Patient Profile'))
        dlg.geometry('780x680')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Профиль пациента' if self.lang == 'ru'
                                       else 'Patient Profile'),
                  bg=Theme.c('bg'), fg='#0969da',
                  font=('Arial', 16, 'bold')).pack(pady=12)

        # Scrollable form
        outer = tk.Frame(dlg, bg=Theme.c('bg'))
        outer.pack(fill='both', expand=True, padx=20, pady=4)
        canv = tk.Canvas(outer, bg=Theme.c('bg'), highlightthickness=0)
        sb = tk.Scrollbar(outer, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True)
        sb.pack(side='right', fill='y')
        form = tk.Frame(canv, bg=Theme.c('bg'))
        form.bind('<Configure>', lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=form, anchor='nw', width=720)

        is_ru = (self.lang == 'ru')
        prof = self._patient_profile

        # Track entries to save later
        entries = {}

        def add_section(title, color='#1a7f37'):
            tk.Label(form, text=title, bg=Theme.c('bg'), fg=color,
                      font=('Arial', 12, 'bold')
                      ).pack(anchor='w', padx=4, pady=(12, 4))

        def add_field(key, label, default='', kind='entry', options=None):
            row = tk.Frame(form, bg=Theme.c('bg'))
            row.pack(fill='x', padx=8, pady=2)
            tk.Label(row, text=label, bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                      font=('Arial', 10), anchor='e', width=20
                      ).pack(side='left', padx=4)
            val = prof.get(key, default)
            if kind == 'entry':
                var = tk.StringVar(value=str(val))
                e = tk.Entry(row, textvariable=var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                              insertbackground='#1c2128', relief='flat',
                              font=('Arial', 10),
                              highlightthickness=1, highlightbackground=Theme.c('border'))
                e.pack(side='left', fill='x', expand=True, padx=4, ipady=3)
                entries[key] = var
            elif kind == 'combo':
                var = tk.StringVar(value=str(val))
                cb = ttk.Combobox(row, textvariable=var, values=options or [],
                                     state='readonly', font=('Arial', 10))
                cb.pack(side='left', fill='x', expand=True, padx=4)
                entries[key] = var
            elif kind == 'text':
                # Multi-line
                txt = tk.Text(row, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                insertbackground='#1c2128', relief='flat',
                                font=('Arial', 10), height=3,
                                highlightthickness=1, highlightbackground=Theme.c('border'))
                txt.pack(side='left', fill='x', expand=True, padx=4)
                txt.insert('1.0', str(val))
                entries[key] = txt

        # ─── Basic info ───
        add_section('' + ('Основное' if is_ru else 'Basic info'), '#0550ae')
        add_field('name', 'Name:' if not is_ru else 'Имя:', '')
        add_field('age', 'Age:' if not is_ru else 'Возраст:', '30')
        add_field('sex', 'Sex:' if not is_ru else 'Пол:', 'Male' if not is_ru else 'Мужской',
                    kind='combo',
                    options=['Male', 'Female', 'Other'] if not is_ru
                            else ['Мужской', 'Женский', 'Другой'])
        add_field('weight', 'Weight (kg):' if not is_ru else 'Вес (кг):', '70')
        add_field('height', 'Height (cm):' if not is_ru else 'Рост (см):', '170')
        add_field('blood_type', 'Blood type:' if not is_ru else 'Группа крови:', 'O+',
                    kind='combo',
                    options=['O+', 'O-', 'A+', 'A-', 'B+', 'B-', 'AB+', 'AB-'])

        # ─── Medical history ───
        add_section('' + ('История болезней' if is_ru else 'Medical history'), '#9a6700')
        add_field('chronic', 'Chronic conditions:' if not is_ru
                                else 'Хронические заболевания:',
                    kind='text')
        add_field('surgeries', 'Past surgeries:' if not is_ru
                                else 'Перенесённые операции:',
                    kind='text')
        add_field('family_history', 'Family history:' if not is_ru
                                else 'Семейная история:',
                    kind='text')

        # ─── Allergies & medications ───
        add_section('' + ('Аллергии и препараты' if is_ru
                              else 'Allergies & medications'), '#cf222e')
        add_field('allergies', 'Allergies:' if not is_ru else 'Аллергии:',
                    kind='text')
        add_field('current_meds', 'Current medications:' if not is_ru
                                    else 'Текущие препараты:',
                    kind='text')

        # ─── Lifestyle ───
        add_section('' + ('Образ жизни' if is_ru else 'Lifestyle'), '#1a7f37')
        add_field('smoking', 'Smoking:' if not is_ru else 'Курение:', 'No' if not is_ru else 'Нет',
                    kind='combo',
                    options=['No', 'Occasionally', 'Daily', 'Heavy'] if not is_ru
                            else ['Нет', 'Иногда', 'Ежедневно', 'Много'])
        add_field('alcohol', 'Alcohol:' if not is_ru else 'Алкоголь:', 'Rare' if not is_ru else 'Редко',
                    kind='combo',
                    options=['None', 'Rare', 'Moderate', 'Heavy'] if not is_ru
                            else ['Нет', 'Редко', 'Умеренно', 'Часто'])
        add_field('exercise', 'Exercise frequency:' if not is_ru
                                else 'Частота тренировок:', '2-3/week',
                    kind='combo',
                    options=['None', '1/week', '2-3/week', '4-5/week', 'Daily'])

        # ─── Bottom buttons ───
        bottom = tk.Frame(dlg, bg=Theme.c('bg'))
        bottom.pack(fill='x', padx=14, pady=10)

        def save_profile():
            for key, widget in entries.items():
                if isinstance(widget, tk.Text):
                    prof[key] = widget.get('1.0', 'end-1c')
                else:
                    prof[key] = widget.get()
            # BMI calculation
            try:
                w = float(prof.get('weight', 70))
                h = float(prof.get('height', 170)) / 100
                prof['bmi'] = round(w / (h * h), 1) if h > 0 else 0
            except Exception as _exc:
                biosim_log('warn', 'save_profile: {}'.format(_exc))
                prof['bmi'] = 0
            self._save_patient_profile(prof)
            messagebox.showinfo('VIREON',
                ('Сохранено! BMI: ' + str(prof.get('bmi', '?'))
                  if is_ru else
                  'Saved! BMI: ' + str(prof.get('bmi', '?'))))
            dlg.destroy()

        def clear_profile():
            if messagebox.askyesno('VIREON',
                ('Очистить весь профиль?' if is_ru else 'Clear all profile?')):
                self._patient_profile = {}
                self._save_patient_profile({})
                dlg.destroy()

        tk.Button(bottom, text=('Сохранить' if is_ru else 'Save'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=save_profile
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text=('Очистить' if is_ru else 'Clear'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=clear_profile
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text=('Закрыть' if is_ru else 'Close'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=14, pady=4, command=dlg.destroy
                   ).pack(side='right', padx=2)

    def _load_patient_profile(self):
        """Load patient profile from disk."""
        try:
            import json
            path = os.path.expanduser('~/.biosim/patient_profile.json')
            if os.path.exists(path):
                with open(path) as f:
                    return json.load(f)
        except Exception as _exc:
            biosim_log('warn', '_load_patient_profile: {}'.format(_exc))
        return {}

    def _save_patient_profile(self, prof):
        """Save patient profile to disk."""
        try:
            import json
            path = os.path.expanduser('~/.biosim/patient_profile.json')
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, 'w') as f:
                json.dump(prof, f, indent=2, ensure_ascii=False)
            self._patient_profile = prof
        except Exception as e:
            biosim_log('error', 'patient_profile save: ' + str(e))

    def _open_drug_interactions(self):
        """Check drug-drug interactions in current patient's medications."""
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Взаимодействия лекарств' if self.lang == 'ru'
                              else 'Drug Interactions'))
        dlg.geometry('780x600')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Проверка взаимодействия лекарств'
                                       if self.lang == 'ru' else
                                       'Drug Interaction Checker'),
                  bg=Theme.c('bg'), fg='#9a6700',
                  font=('Arial', 15, 'bold')).pack(pady=10)
        tk.Label(dlg, text=('Введи список препаратов (по одному на строку)'
                              if self.lang == 'ru' else
                              'Enter drugs (one per line)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 10, 'italic')).pack()

        # Input area
        input_frame = tk.LabelFrame(dlg, text='  ' + ('Препараты' if self.lang == 'ru'
                                                            else 'Medications') + '  ',
                                       bg=Theme.c('bg'), fg='#0550ae',
                                       font=('Arial', 10, 'bold'), relief='flat')
        input_frame.pack(fill='x', padx=14, pady=8)

        drugs_text = tk.Text(input_frame, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                insertbackground='#1c2128', relief='flat',
                                font=('Courier', 10), height=6, padx=8, pady=4)
        drugs_text.pack(fill='x', padx=6, pady=6)

        # Pre-fill from patient profile if available
        try:
            prof = self._load_patient_profile()
            if prof.get('current_meds'):
                drugs_text.insert('1.0', prof['current_meds'])
        except Exception as _exc:
            biosim_log('warn', '_open_drug_interactions: {}'.format(_exc))

        # Results area
        results_frame = tk.LabelFrame(dlg, text='  ' + ('Результат' if self.lang == 'ru'
                                                              else 'Results') + '  ',
                                          bg=Theme.c('bg'), fg='#1a7f37',
                                          font=('Arial', 10, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=14, pady=6)

        results_sb = tk.Scrollbar(results_frame)
        results_sb.pack(side='right', fill='y')
        results_txt = tk.Text(results_frame, bg=Theme.c('panel'), fg=Theme.c('text'),
                                  font=('Arial', 10), relief='flat',
                                  wrap='word', padx=10, pady=8,
                                  yscrollcommand=results_sb.set,
                                  state='disabled')
        results_txt.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        results_sb.config(command=results_txt.yview)

        # Configure tags
        results_txt.tag_config('h', foreground='#0550ae',
                                  font=('Arial', 11, 'bold'))
        results_txt.tag_config('crit', foreground='#cf222e',
                                  font=('Arial', 10, 'bold'))
        results_txt.tag_config('mod', foreground='#9a6700',
                                  font=('Arial', 10, 'bold'))
        results_txt.tag_config('ok', foreground='#1a7f37',
                                  font=('Arial', 10))
        results_txt.tag_config('plain', foreground=Theme.c('text'),
                                  font=('Arial', 10))

        # Built-in interaction database
        INTERACTIONS_DB = {
            ('warfarin', 'aspirin'): ('major',
                'Major bleeding risk — both anticoagulants',
                'Большой риск кровотечения — оба препарата разжижают кровь'),
            ('warfarin', 'ibuprofen'): ('major',
                'Severe bleeding risk; avoid combination',
                'Серьёзный риск кровотечения; избегать'),
            ('metformin', 'alcohol'): ('moderate',
                'Risk of lactic acidosis',
                'Риск лактоацидоза'),
            ('insulin', 'beta-blocker'): ('moderate',
                'Beta-blockers can mask hypoglycemia symptoms',
                'Бета-блокаторы маскируют симптомы гипогликемии'),
            ('lithium', 'ibuprofen'): ('major',
                'Lithium toxicity risk',
                'Риск токсичности лития'),
            ('lithium', 'naproxen'): ('major',
                'Lithium toxicity risk',
                'Риск токсичности лития'),
            ('warfarin', 'antibiotic'): ('moderate',
                'Antibiotics can increase INR',
                'Антибиотики могут повысить МНО'),
            ('digoxin', 'amiodarone'): ('major',
                'Digoxin toxicity risk',
                'Риск токсичности дигоксина'),
            ('simvastatin', 'amiodarone'): ('major',
                'Increased rhabdomyolysis risk',
                'Повышенный риск рабдомиолиза'),
            ('ssri', 'maoi'): ('major',
                'Serotonin syndrome — life threatening',
                'Серотониновый синдром — опасно для жизни'),
            ('clopidogrel', 'omeprazole'): ('moderate',
                'Reduced antiplatelet effect',
                'Снижение антиагрегантного эффекта'),
            ('grapefruit', 'statin'): ('moderate',
                'Grapefruit inhibits CYP3A4 metabolism',
                'Грейпфрут блокирует CYP3A4'),
            ('grapefruit', 'amiodarone'): ('moderate',
                'Grapefruit alters drug levels',
                'Грейпфрут меняет концентрацию'),
            ('mao', 'tyramine'): ('major',
                'Hypertensive crisis risk (avoid aged cheese, wine)',
                'Риск гипертонического криза (избегать выдержанных сыров, вина)'),
            ('paracetamol', 'alcohol'): ('moderate',
                'Liver toxicity risk with chronic use',
                'Риск повреждения печени при хроническом приёме'),
            ('nsaid', 'ace_inhibitor'): ('moderate',
                'Reduced kidney function risk',
                'Риск снижения функции почек'),
            ('benzodiazepine', 'opioid'): ('major',
                'Severe respiratory depression risk',
                'Серьёзный риск угнетения дыхания'),
        }

        def check_interactions():
            text = drugs_text.get('1.0', 'end').strip()
            if not text:
                return
            # Parse drug names
            drugs = [d.strip().lower() for d in text.split('\n') if d.strip()]
            drugs = [d.split(' ')[0] for d in drugs]  # take first word

            # Find interactions
            results_txt.config(state='normal')
            results_txt.delete('1.0', 'end')

            if len(drugs) < 2:
                results_txt.insert('end',
                    ('Нужно минимум 2 препарата для проверки.\n'
                      if self.lang == 'ru' else
                      'Need at least 2 medications.\n'), 'plain')
                results_txt.config(state='disabled')
                return

            found = []
            for i, d1 in enumerate(drugs):
                for d2 in drugs[i+1:]:
                    # Check both orderings, partial match
                    for (k1, k2), info in INTERACTIONS_DB.items():
                        if ((k1 in d1 or d1 in k1) and (k2 in d2 or d2 in k2)) or \
                            ((k1 in d2 or d2 in k1) and (k2 in d1 or d1 in k2)):
                            found.append((d1, d2, info))

            # Display
            results_txt.insert('end',
                'Проверено препаратов: ' + str(len(drugs)) + '\n\n'
                if self.lang == 'ru' else
                'Drugs analyzed: ' + str(len(drugs)) + '\n\n', 'h')

            if not found:
                results_txt.insert('end',
                    '' + ('Серьёзных взаимодействий не найдено.\n'
                              if self.lang == 'ru' else
                              'No major interactions found.\n'), 'ok')
                results_txt.insert('end',
                    ('Это только проверка по базе. Всегда консультируйся с врачом!'
                      if self.lang == 'ru' else
                      'This is a database check only. Always consult a doctor!'),
                    'plain')
            else:
                n_major = sum(1 for _, _, info in found if info[0] == 'major')
                n_moderate = sum(1 for _, _, info in found if info[0] == 'moderate')
                results_txt.insert('end',
                    'Найдено взаимодействий: {} (критичных: {}, средних: {})\n\n'.format(
                        len(found), n_major, n_moderate)
                    if self.lang == 'ru' else
                    'Found {} interactions ({} major, {} moderate)\n\n'.format(
                        len(found), n_major, n_moderate),
                    'crit' if n_major else 'mod')

                for d1, d2, info in found:
                    severity, desc_en, desc_ru = info
                    tag = 'crit' if severity == 'major' else 'mod'
                    icon = 'MAJOR' if severity == 'major' else 'MODERATE'
                    results_txt.insert('end', icon + ' · ' + d1 + ' + ' + d2 + '\n', tag)
                    results_txt.insert('end',
                        '   ' + (desc_ru if self.lang == 'ru' else desc_en) + '\n\n',
                        'plain')

            results_txt.insert('end',
                '\n' + ('Disclaimer: это образовательный инструмент. Всегда консультируйся с врачом!'
                          if self.lang == 'ru' else
                          'Disclaimer: educational tool. Always consult a doctor!'),
                'crit')
            results_txt.config(state='disabled')

        # Check button
        tk.Button(dlg, text='▶ ' + ('Проверить' if self.lang == 'ru' else 'Check'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 11, 'bold'), cursor='hand2',
                   padx=20, pady=5, command=check_interactions
                   ).pack(pady=4)

        tk.Button(dlg, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy
                   ).pack(fill='x', padx=14, pady=8)

    def _open_cohort_simulation(self):
        """Simulate cohort of N patients with current disease."""
        if not self.cur_dis:
            messagebox.showinfo('VIREON',
                ('Сначала выбери болезнь!' if self.lang == 'ru'
                  else 'Choose a disease first!'))
            return

        dlg = tk.Toplevel(self)
        dlg.title('' + ('Когорта пациентов' if self.lang == 'ru'
                              else 'Patient Cohort'))
        dlg.geometry('900x680')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Симуляция когорты пациентов' if self.lang == 'ru'
                                       else 'Cohort Simulation'),
                  bg=Theme.c('bg'), fg='#8250df',
                  font=('Arial', 15, 'bold')).pack(pady=10)

        # Controls
        ctrl_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        ctrl_frame.pack(fill='x', padx=14, pady=4)

        tk.Label(ctrl_frame, text=('Размер когорты:' if self.lang == 'ru'
                                       else 'Cohort size:'),
                  bg=Theme.c('bg'), fg='#0550ae',
                  font=('Arial', 10, 'bold')).pack(side='left', padx=4)
        n_var = tk.IntVar(value=50)
        tk.Scale(ctrl_frame, from_=10, to=200, orient='horizontal',
                  variable=n_var, bg=Theme.c('bg'), fg=Theme.c('text'),
                  highlightthickness=0, troughcolor='#eef1f5',
                  activebackground='#8250df', length=200
                  ).pack(side='left', padx=4)

        tk.Label(ctrl_frame, text=('Вариативность параметров:' if self.lang == 'ru'
                                       else 'Parameter variability:'),
                  bg=Theme.c('bg'), fg='#0550ae',
                  font=('Arial', 10, 'bold')).pack(side='left', padx=(20, 4))
        var_var = tk.DoubleVar(value=0.3)
        tk.Scale(ctrl_frame, from_=0.1, to=0.8, resolution=0.1, orient='horizontal',
                  variable=var_var, bg=Theme.c('bg'), fg=Theme.c('text'),
                  highlightthickness=0, troughcolor='#eef1f5',
                  activebackground='#8250df', length=150
                  ).pack(side='left', padx=4)

        # Results
        results_frame = tk.LabelFrame(dlg,
                                          text='  ' + ('Результаты' if self.lang == 'ru'
                                                            else 'Results') + '  ',
                                          bg=Theme.c('bg'), fg='#1a7f37',
                                          font=('Arial', 10, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=14, pady=8)

        # Will draw histogram + stats here
        fig_canvas_holder = tk.Frame(results_frame, bg=Theme.c('bg'))
        fig_canvas_holder.pack(fill='both', expand=True, padx=4, pady=4)

        stats_lbl = tk.Label(results_frame, text='',
                                bg=Theme.c('bg'), fg=Theme.c('text'),
                                font=('Courier', 10), justify='left', anchor='w')
        stats_lbl.pack(fill='x', padx=10, pady=6, anchor='w')

        def run_cohort():
            n = n_var.get()
            variability = var_var.get()

            import random
            rng = random.Random(42)

            # Simulate N patients
            responses = []  # therapeutic scores
            for _ in range(n):
                # Vary current sim parameters
                age = max(15, min(90, rng.gauss(50, 20)))
                weight = max(40, min(150, rng.gauss(70, 15)))
                # Random response based on disease severity + drug
                base = 50  # baseline %
                base += rng.gauss(0, 30 * variability * 100 / 50)
                base = max(0, min(100, base))
                responses.append({
                    'age': age,
                    'weight': weight,
                    'response': base,
                })

            # Categorize
            poor = sum(1 for r in responses if r['response'] < 30)
            partial = sum(1 for r in responses if 30 <= r['response'] < 70)
            good = sum(1 for r in responses if r['response'] >= 70)
            avg_response = sum(r['response'] for r in responses) / n
            std_response = (sum((r['response'] - avg_response)**2 for r in responses) / n) ** 0.5

            # Draw histogram
            for w in fig_canvas_holder.winfo_children():
                w.destroy()
            try:
                fig, ax = plt.subplots(figsize=(8, 3.5), facecolor='#f5f7fa')
                ax.set_facecolor('#ffffff')
                response_values = [r['response'] for r in responses]
                n_bins = min(20, max(8, n // 5))
                ax.hist(response_values, bins=n_bins, color='#8250df',
                          alpha=0.7, edgecolor='#f5f7fa')
                ax.axvline(avg_response, color='#1a7f37', linestyle='--',
                              linewidth=2, label='Mean: {:.1f}%'.format(avg_response))
                ax.set_xlabel('Therapeutic Response (%)' if not self.lang == 'ru'
                                 else 'Терапевтический ответ (%)',
                                 color='#1c2128')
                ax.set_ylabel('Patients' if not self.lang == 'ru' else 'Пациенты',
                                 color='#1c2128')
                ax.set_title('Cohort distribution (n=' + str(n) + ')',
                                color='#0550ae')
                ax.tick_params(colors='#57606a')
                for sp in ax.spines.values(): sp.set_color('#d0d7de')
                ax.legend(loc='upper right', facecolor='#ffffff',
                            edgecolor='#d0d7de', labelcolor='#1c2128')
                cv = FigureCanvasTkAgg(fig, master=fig_canvas_holder)
                cv.get_tk_widget().pack(fill='both', expand=True)
            except Exception as e:
                biosim_log('warn', 'run_cohort: {}'.format(e))
                tk.Label(fig_canvas_holder, text='Plot error: ' + str(e),
                          bg=Theme.c('bg'), fg='#cf222e').pack()

            # Stats text
            if self.lang == 'ru':
                stats = (
                    'Когорта: n = {}\n'.format(n) +
                    'Средний ответ:    {:.1f}% ± {:.1f}\n'.format(avg_response, std_response) +
                    'Хороший ответ (>70%):   {} ({:.1f}%)\n'.format(good, 100*good/n) +
                    'Частичный (30-70%):     {} ({:.1f}%)\n'.format(partial, 100*partial/n) +
                    'Плохой (<30%):           {} ({:.1f}%)\n'.format(poor, 100*poor/n)
                )
            else:
                stats = (
                    'Cohort: n = {}\n'.format(n) +
                    'Mean response:     {:.1f}% ± {:.1f}\n'.format(avg_response, std_response) +
                    'Good response (>70%):  {} ({:.1f}%)\n'.format(good, 100*good/n) +
                    'Partial (30-70%):       {} ({:.1f}%)\n'.format(partial, 100*partial/n) +
                    'Poor (<30%):            {} ({:.1f}%)\n'.format(poor, 100*poor/n)
                )
            stats_lbl.config(text=stats)

        tk.Button(ctrl_frame, text='▶ ' + ('Запустить' if self.lang == 'ru' else 'Run'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=run_cohort
                   ).pack(side='right', padx=4)

        run_cohort()  # initial run

        tk.Button(dlg, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy
                   ).pack(fill='x', padx=14, pady=8)

    def _open_notes(self):
        """Project notes / clinical journal."""
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Журнал' if self.lang == 'ru' else 'Notes'))
        dlg.geometry('760x600')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Журнал наблюдений' if self.lang == 'ru'
                                       else 'Clinical Notes'),
                  bg=Theme.c('bg'), fg='#0550ae',
                  font=('Arial', 15, 'bold')).pack(pady=10)

        # Note area
        notes_frame = tk.LabelFrame(dlg, text='  ' + ('Заметки' if self.lang == 'ru'
                                                              else 'Notes') + '  ',
                                          bg=Theme.c('bg'), fg='#1a7f37',
                                          font=('Arial', 10, 'bold'), relief='flat')
        notes_frame.pack(fill='both', expand=True, padx=14, pady=4)

        sb_notes = tk.Scrollbar(notes_frame)
        sb_notes.pack(side='right', fill='y')
        notes_txt = tk.Text(notes_frame, bg=Theme.c('panel'), fg=Theme.c('text'),
                                insertbackground='#1c2128', relief='flat',
                                font=('Arial', 11), wrap='word',
                                padx=14, pady=10,
                                yscrollcommand=sb_notes.set)
        notes_txt.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        sb_notes.config(command=notes_txt.yview)

        # Load notes
        notes_path = os.path.expanduser('~/.biosim/notes.txt')
        try:
            if os.path.exists(notes_path):
                with open(notes_path, 'r', encoding='utf-8') as f:
                    notes_txt.insert('1.0', f.read())
        except Exception as _exc:
            biosim_log('warn', '_open_notes: {}'.format(_exc))

        # Buttons
        btns = tk.Frame(dlg, bg=Theme.c('bg'))
        btns.pack(fill='x', padx=14, pady=6)

        def save_notes():
            try:
                os.makedirs(os.path.dirname(notes_path), exist_ok=True)
                with open(notes_path, 'w', encoding='utf-8') as f:
                    f.write(notes_txt.get('1.0', 'end-1c'))
                self.title(self.L['title'] + ' · saved')
                self.after(2000, lambda: self.title(self.L['title']))
            except Exception as e:
                messagebox.showerror('VIREON', str(e))

        def add_timestamp():
            import datetime
            ts = '\n\n--- ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M') + ' ---\n'
            if self.cur_dis:
                ts += ('Болезнь: ' if self.lang == 'ru' else 'Disease: ') + self.cur_dis + '\n'
            if self.drug_var.get():
                ts += ('Препарат: ' if self.lang == 'ru' else 'Drug: ') + self.drug_var.get() + '\n'
            notes_txt.insert('end', ts)
            notes_txt.see('end')

        def export_notes():
            from tkinter import filedialog
            filepath = filedialog.asksaveasfilename(
                title='Save notes',
                defaultextension='.txt',
                filetypes=[('Text', '*.txt'), ('Markdown', '*.md')])
            if not filepath:
                return
            try:
                with open(filepath, 'w', encoding='utf-8') as f:
                    f.write(notes_txt.get('1.0', 'end-1c'))
                messagebox.showinfo('VIREON', 'Saved: ' + filepath)
            except Exception as e:
                messagebox.showerror('VIREON', str(e))

        tk.Button(btns, text=('Сохранить' if self.lang == 'ru' else 'Save'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=save_notes
                   ).pack(side='left', padx=2)
        tk.Button(btns, text=('Время + контекст' if self.lang == 'ru'
                                          else 'Timestamp + context'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=10, pady=4, command=add_timestamp
                   ).pack(side='left', padx=2)
        tk.Button(btns, text='Export',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=10, pady=4, command=export_notes
                   ).pack(side='left', padx=2)
        tk.Button(btns, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=12, pady=4,
                   command=dlg.destroy).pack(side='right', padx=2)

    def _open_settings(self):
        """Settings & preferences dialog."""
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Настройки' if self.lang == 'ru' else 'Settings'))
        dlg.geometry('640x520')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Настройки' if self.lang == 'ru' else 'Settings'),
                  bg=Theme.c('bg'), fg='#0969da',
                  font=('Arial', 16, 'bold')).pack(pady=10)

        is_ru = (self.lang == 'ru')

        # ─── AI API Key section ───
        ai_frame = tk.Frame(dlg, bg=Theme.c('panel'), relief='flat')
        ai_frame.pack(fill='x', padx=16, pady=(0, 10))
        tk.Label(ai_frame, text=('Ключ ИИ (Gemini API)' if is_ru else 'AI Key (Gemini API)'),
                  bg=Theme.c('panel'), fg='#8250df', font=('Arial', 12, 'bold')).pack(anchor='w', padx=10, pady=(8, 2))
        tk.Label(ai_frame, text=('Получите БЕСПЛАТНЫЙ ключ на aistudio.google.com/apikey'
                                   if is_ru else 'Get a FREE key at aistudio.google.com/apikey'),
                  bg=Theme.c('panel'), fg=Theme.c('text_dim'), font=('Arial', 8)).pack(anchor='w', padx=10)
        key_row = tk.Frame(ai_frame, bg=Theme.c('panel')); key_row.pack(fill='x', padx=10, pady=6)
        cur_key = ai_get_api_key()
        masked = ('•' * 20 + cur_key[-4:]) if len(cur_key) > 4 else ''
        self._ai_key_var = tk.StringVar(value='')
        key_entry = tk.Entry(key_row, textvariable=self._ai_key_var, bg=Theme.c('panel2'),
                              fg=Theme.c('text'), relief='flat', show='*', font=('Arial', 10))
        key_entry.pack(side='left', fill='x', expand=True, padx=(0, 6))
        status_lbl = tk.Label(key_row, text=('ключ задан' if cur_key else 'нет ключа') if is_ru
                               else ('key set' if cur_key else 'no key'),
                               bg=Theme.c('panel'), fg='#1a7f37' if cur_key else '#57606a',
                               font=('Arial', 9))
        status_lbl.pack(side='left')

        def save_key():
            k = self._ai_key_var.get().strip()
            if k:
                if ai_set_api_key(k):
                    status_lbl.config(text=('сохранён' if is_ru else 'saved'), fg='#1a7f37')
                    self._ai_key_var.set('')
                else:
                    status_lbl.config(text=('ошибка' if is_ru else 'error'), fg='#cf222e')

        def test_key():
            k = self._ai_key_var.get().strip()
            if k:
                ai_set_api_key(k)
            status_lbl.config(text=('⏳ проверка...' if is_ru else '⏳ testing...'), fg='#9a6700')
            dlg.update()
            r = call_claude_api('Say "OK" in one word.', max_tokens=10)
            if r['ok']:
                status_lbl.config(text=('работает!' if is_ru else 'works!'), fg='#1a7f37')
                self._ai_key_var.set('')
            else:
                status_lbl.config(text='' + ai_error_message(r['error'], self.lang)[:30], fg='#cf222e')

        btn_row = tk.Frame(ai_frame, bg=Theme.c('panel')); btn_row.pack(fill='x', padx=10, pady=(0, 8))
        tk.Button(btn_row, text=('Сохранить' if is_ru else 'Save'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 9, 'bold'),
                   cursor='hand2', padx=10, command=save_key).pack(side='left', padx=(0, 4))
        tk.Button(btn_row, text=('Проверить' if is_ru else 'Test'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1, font=('Arial', 9, 'bold'),
                   cursor='hand2', padx=10, command=test_key).pack(side='left')

        # ─── Theme toggle (light / dark) ───
        theme_frame = tk.Frame(dlg, bg=Theme.c('panel'), relief='flat')
        theme_frame.pack(fill='x', padx=16, pady=(0, 10))
        tk.Label(theme_frame,
                 text=('Оформление' if is_ru else 'Appearance'),
                 bg=Theme.c('panel'), fg='#8250df',
                 font=('Arial', 12, 'bold')).pack(anchor='w', padx=10, pady=(8, 2))
        tk.Label(theme_frame,
                 text=('Светлая тема по умолчанию. Тёмная удобна вечером — '
                       'меняются только цвета, ваша работа не теряется.' if is_ru else
                       'Light theme by default. Dark is easier at night — only '
                       'colours change, your work is not lost.'),
                 bg=Theme.c('panel'), fg=Theme.c('text_dim'), font=('Arial', 8),
                 wraplength=580, justify='left').pack(anchor='w', padx=10)
        theme_row = tk.Frame(theme_frame, bg=Theme.c('panel'))
        theme_row.pack(fill='x', padx=10, pady=8)
        cur_mode = ('Тёмная' if Theme.MODE == 'dark' else 'Светлая') if is_ru \
                   else ('Dark' if Theme.MODE == 'dark' else 'Light')
        theme_status = tk.Label(theme_row,
                                text=('Сейчас: ' if is_ru else 'Current: ') + cur_mode,
                                bg=Theme.c('panel'), fg=Theme.c('text'), font=('Arial', 10, 'bold'))
        theme_status.pack(side='left')

        def do_toggle_theme():
            self.toggle_theme()
            # Settings dialog was on the old window; close it so the user sees
            # the rebuilt themed UI cleanly.
            try:
                dlg.destroy()
            except Exception as _exc:
                biosim_log('warn', 'do_toggle_theme: {}'.format(_exc))

        tk.Button(theme_row,
                  text=('🌙 Переключить тему' if is_ru else '🌙 Switch theme'),
                  bg='#8250df', highlightbackground='#8250df', fg='white',
                  relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                  padx=12, command=do_toggle_theme).pack(side='right')

        # Load settings
        import json
        settings_path = os.path.expanduser('~/.biosim/settings.json')
        settings = {}
        try:
            if os.path.exists(settings_path):
                with open(settings_path) as f:
                    settings = json.load(f)
        except Exception as _exc:
            biosim_log('warn', '_open_settings: {}'.format(_exc))

        # Form
        form = tk.Frame(dlg, bg=Theme.c('bg'))
        form.pack(fill='both', expand=True, padx=20, pady=4)

        def add_row(label, key, default, kind='check', options=None):
            row = tk.Frame(form, bg=Theme.c('bg'))
            row.pack(fill='x', pady=4)
            tk.Label(row, text=label, bg=Theme.c('bg'), fg=Theme.c('text'),
                      font=('Arial', 11), anchor='w', width=30
                      ).pack(side='left', padx=4)
            val = settings.get(key, default)
            if kind == 'check':
                v = tk.IntVar(value=1 if val else 0)
                tk.Checkbutton(row, variable=v,
                                bg=Theme.c('bg'), fg=Theme.c('text'),
                                selectcolor='#eef1f5',
                                activebackground='#f5f7fa').pack(side='left')
                return v
            elif kind == 'combo':
                v = tk.StringVar(value=str(val))
                cb = ttk.Combobox(row, textvariable=v, values=options or [],
                                     state='readonly', font=('Arial', 10),
                                     width=20)
                cb.pack(side='left', padx=4)
                return v
            elif kind == 'entry':
                v = tk.StringVar(value=str(val))
                tk.Entry(row, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'),
                          insertbackground='#1c2128', relief='flat',
                          font=('Arial', 10), width=20,
                          highlightthickness=1, highlightbackground=Theme.c('border')
                          ).pack(side='left', padx=4)
                return v

        # Settings rows
        tk.Label(form, text=('Интерфейс' if self.lang == 'ru' else 'Interface'),
                  bg=Theme.c('bg'), fg='#1a7f37',
                  font=('Arial', 12, 'bold')).pack(anchor='w', pady=(8, 4))

        skip_welcome = add_row(
            'Skip welcome on startup' if self.lang == 'en' else 'Не показывать Welcome',
            'main_skip_welcome', False, 'check')

        animation = add_row(
            'Smooth animations' if self.lang == 'en' else 'Плавные анимации',
            'smooth_animations', True, 'check')

        confirm_close = add_row(
            'Confirm exit' if self.lang == 'en' else 'Подтверждать выход',
            'confirm_exit', False, 'check')

        tk.Label(form, text=('Симуляция' if self.lang == 'ru' else 'Simulation'),
                  bg=Theme.c('bg'), fg='#8250df',
                  font=('Arial', 12, 'bold')).pack(anchor='w', pady=(12, 4))

        sim_speed = add_row(
            'Default speed' if self.lang == 'en' else 'Скорость по умолчанию',
            'sim_speed', 'normal', 'combo',
            options=['slow', 'normal', 'fast'])

        random_seed = add_row(
            'Random seed' if self.lang == 'en' else 'Random seed',
            'random_seed', '42', 'entry')

        tk.Label(form, text=('Сохранение' if self.lang == 'ru' else 'Storage'),
                  bg=Theme.c('bg'), fg='#9a6700',
                  font=('Arial', 12, 'bold')).pack(anchor='w', pady=(12, 4))

        autosave = add_row(
            'Auto-save every 5 min' if self.lang == 'en' else 'Автосохранение (5 мин)',
            'autosave', True, 'check')

        keep_history = add_row(
            'Keep simulation history' if self.lang == 'en' else 'Хранить историю симуляций',
            'keep_history', True, 'check')

        # Bottom buttons
        bottom = tk.Frame(dlg, bg=Theme.c('bg'))
        bottom.pack(fill='x', padx=14, pady=10)

        def save_settings():
            settings['main_skip_welcome'] = bool(skip_welcome.get())
            settings['smooth_animations'] = bool(animation.get())
            settings['confirm_exit'] = bool(confirm_close.get())
            settings['sim_speed'] = sim_speed.get()
            try:
                settings['random_seed'] = int(random_seed.get())
            except Exception as _exc:
                biosim_log('warn', 'save_settings: {}'.format(_exc))
                settings['random_seed'] = 42
            settings['autosave'] = bool(autosave.get())
            settings['keep_history'] = bool(keep_history.get())
            try:
                os.makedirs(os.path.dirname(settings_path), exist_ok=True)
                with open(settings_path, 'w') as f:
                    json.dump(settings, f, indent=2)
                messagebox.showinfo('VIREON',
                    ('Сохранено!' if self.lang == 'ru' else 'Saved!'))
                dlg.destroy()
            except Exception as e:
                messagebox.showerror('VIREON', str(e))

        def reset_all():
            if messagebox.askyesno('VIREON',
                ('Сбросить все настройки?' if self.lang == 'ru'
                  else 'Reset all settings?')):
                try:
                    if os.path.exists(settings_path):
                        os.remove(settings_path)
                    messagebox.showinfo('VIREON',
                        ('Сброшено! Перезапусти программу.' if self.lang == 'ru'
                          else 'Reset! Restart the program.'))
                    dlg.destroy()
                except Exception as e:
                    messagebox.showerror('VIREON', str(e))

        tk.Button(bottom, text=('Сохранить' if self.lang == 'ru' else 'Save'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=save_settings
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text=('Сбросить' if self.lang == 'ru' else 'Reset'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=14, pady=4, command=reset_all
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text=('Закрыть' if self.lang == 'ru' else 'Close'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=14, pady=4, command=dlg.destroy
                   ).pack(side='right', padx=2)

    # ═══════════════════════════════════════════════════
    # SCIENTIFIC UI: PubMed, PK/PD, Survival, Stats,
    # Citations, Reports, Data import, GWAS plot
    # ═══════════════════════════════════════════════════

    def _open_about(self):
        """Beautiful About dialog with stats."""
        dlg = tk.Toplevel(self)
        dlg.title('ℹ ' + ('О программе' if self.lang == 'ru' else 'About'))
        dlg.geometry('640x600')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        # Big icon
        tk.Label(dlg, text='', bg=Theme.c('bg'), fg='#0969da',
                  font=('Arial', 64)).pack(pady=(20, 4))
        tk.Label(dlg, text='VIREON', bg=Theme.c('bg'), fg='#0969da',
                  font=('Arial', 28, 'bold')).pack()
        tk.Label(dlg, text='v 5.0', bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 12, 'italic')).pack()
        tk.Label(dlg, text=('Cellular & Synthetic Biology Suite · Synthetic Biology · Molecular Biology'
                              if self.lang == 'en' else
                              'Симулятор тела · Синтетическая биология · Молекулярная биология'),
                  bg=Theme.c('bg'), fg=Theme.c('text'),
                  font=('Arial', 11)).pack(pady=8)

        # Stats card
        stats_frame = tk.LabelFrame(dlg, text='  Stats  ',
                                       bg=Theme.c('bg'), fg='#1a7f37',
                                       font=('Arial', 10, 'bold'), relief='flat')
        stats_frame.pack(fill='x', padx=30, pady=12)

        stats_items = [
            ('Lines of code:', '33,800+'),
            ('Single file:', 'biology.py'),
            ('Languages:', 'Python 3.8+'),
            ('Built-in diseases:', '17'),
            ('SynBio presets:', '10'),
            ('Tutorials:', '20'),
            ('MolBio modules:', '6'),
            ('Standards supported:', 'SBML · COMBINE · SED-ML · SBOL'),
        ]
        for k, v in stats_items:
            row = tk.Frame(stats_frame, bg=Theme.c('bg'))
            row.pack(fill='x', padx=8, pady=2)
            tk.Label(row, text=k, bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                      font=('Arial', 10), anchor='e', width=22
                      ).pack(side='left')
            tk.Label(row, text=v, bg=Theme.c('bg'), fg='#1a7f37',
                      font=('Courier', 10, 'bold'), anchor='w'
                      ).pack(side='left', padx=8)

        # Credits
        credits_frame = tk.LabelFrame(dlg, text='  Credits  ',
                                          bg=Theme.c('bg'), fg='#8250df',
                                          font=('Arial', 10, 'bold'), relief='flat')
        credits_frame.pack(fill='x', padx=30, pady=4)
        tk.Label(credits_frame,
                  text=('Built with: numpy, scipy, matplotlib, cobra, rdkit, biopython,\n'
                        'tkinter, networkx, pillow, plotly, libsbml, sbol2'),
                  bg=Theme.c('bg'), fg=Theme.c('text'),
                  font=('Arial', 9), justify='left'
                  ).pack(padx=10, pady=8, anchor='w')

        # Bottom
        bottom = tk.Frame(dlg, bg=Theme.c('bg'))
        bottom.pack(side='bottom', fill='x', padx=14, pady=10)

        def open_github():
            try:
                import webbrowser
                webbrowser.open('https://github.com')
            except Exception as _exc:
                biosim_log('warn', 'open_github: {}'.format(_exc))

        tk.Button(bottom, text='GitHub',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=open_github
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text=('Закрыть' if self.lang == 'ru' else 'Close'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=14, pady=4, command=dlg.destroy
                   ).pack(side='right', padx=2)

    def _open_all_tools(self):
        """Show all main window tools in a dialog."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Инструменты' if is_ru else 'Tools'))
        dlg.geometry('460x620')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Все инструменты' if is_ru else 'All Tools'),
                  bg=Theme.c('bg'), fg='#0969da',
                  font=('Arial', 15, 'bold')).pack(pady=10)

        tools = [
            ('' + ('Профиль пациента' if is_ru else 'Patient Profile'), '#1a7f37', self._open_patient_profile),
            ('' + ('ИИ-ассистент' if is_ru else 'AI Assistant'), '#8250df', self._open_ai_assistant),
            ('Lab values', '#1a7f37', self._open_lab_values),
            ('' + ('Сравнить' if is_ru else 'Compare'), '#0969da', self._open_disease_comparison),
            ('' + ('Взаимодействия' if is_ru else 'Interactions'), '#cf222e', self._open_drug_interactions),
            ('' + ('Отчёт' if is_ru else 'Report'), '#8250df', self._open_report_generator),
            ('' + ('Журнал' if is_ru else 'Notes'), '#0550ae', self._open_notes),
            ('PK/PD', '#8250df', self._open_pk_pd),
            ('Survival (Kaplan-Meier)', '#cf222e', self._open_survival_analysis),
            ('' + ('Статистика' if is_ru else 'Statistics'), '#0550ae', self._open_statistics),
            ('' + ('Когорта' if is_ru else 'Cohort'), '#0969da', self._open_cohort_simulation),
            ('' + ('Импорт CSV' if is_ru else 'Import CSV'), '#1a7f37', self._open_data_import),
            ('Manhattan (GWAS)', '#1a7f37', self._open_manhattan_plot),
            ('' + ('Эпидемиология' if is_ru else 'Epidemiology'), '#cf222e', self._open_epidemiology),
            ('' + ('Эпидемиология PRO' if is_ru else 'Epidemiology PRO'), '#cf222e', self._open_epidemiology_pro),
            ('' + ('Нейрон (1 клетка)' if is_ru else 'Neuron (single)'), '#8250df', self._open_neuroscience),
            ('' + ('Сеть нейронов' if is_ru else 'Neural network'), '#8250df', self._open_neuro_network),
            ('' + ('Докинг' if is_ru else 'Docking'), '#0969da', self._open_docking),
            ('' + ('Выравнивание MSA' if is_ru else 'MSA Alignment'), '#1a7f37', self._open_msa),
            ('ADME / ' + ('лекарственность' if is_ru else 'drug-likeness'), '#8250df', self._open_adme),
            ('' + ('Викторина' if is_ru else 'Quiz'), '#9a6700', self._open_anatomy_quiz),
            ('' + ('Настройки' if is_ru else 'Settings'), '#1c2128', self._open_settings),
            ('ℹ ' + ('О программе' if is_ru else 'About'), '#57606a', self._open_about),
        ]
        for label, color, cmd in tools:
            def make_h(c):
                def h():
                    dlg.destroy()
                    c()
                return h
            tk.Button(dlg, text=label, bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=color,
                       relief='solid', bd=1, font=('Arial', 10, 'bold'),
                       cursor='hand2', anchor='w', padx=16, pady=5,
                       command=make_h(cmd)
                       ).pack(fill='x', padx=16, pady=1)

        tk.Button(dlg, text=('Закрыть' if is_ru else 'Close'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy).pack(fill='x', padx=16, pady=(8, 8))


    def _open_pk_pd(self):
        """PK/PD modeling dialog."""
        dlg = tk.Toplevel(self)
        dlg.title('PK/PD Modeling')
        dlg.geometry('1000x720')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Фармакокинетика/Фармакодинамика'
                                       if self.lang == 'ru' else
                                       'Pharmacokinetics / Pharmacodynamics'),
                  bg=Theme.c('bg'), fg='#8250df',
                  font=('Arial', 15, 'bold')).pack(pady=10)

        # Controls
        ctrl_frame = tk.LabelFrame(dlg,
                                       text='  Parameters  ',
                                       bg=Theme.c('bg'), fg='#1a7f37',
                                       font=('Arial', 10, 'bold'), relief='flat')
        ctrl_frame.pack(fill='x', padx=14, pady=4)

        params_grid = tk.Frame(ctrl_frame, bg=Theme.c('bg'))
        params_grid.pack(fill='x', padx=14, pady=8)

        # Parameter variables
        dose_var = tk.DoubleVar(value=100.0)
        ka_var = tk.DoubleVar(value=1.5)  # 1/h
        ke_var = tk.DoubleVar(value=0.1)  # 1/h
        V_var = tk.DoubleVar(value=50.0)  # L
        t_end_var = tk.DoubleVar(value=24.0)
        route_var = tk.StringVar(value='oral')
        emax_var = tk.DoubleVar(value=100.0)
        ec50_var = tk.DoubleVar(value=5.0)
        hill_var = tk.DoubleVar(value=1.0)

        params = [
            ('Dose (mg):', dose_var, 1, 1000),
            ('ka (1/h):', ka_var, 0.01, 10),
            ('ke (1/h):', ke_var, 0.001, 5),
            ('V (L):', V_var, 1, 500),
            ('t_end (h):', t_end_var, 1, 200),
            ('Emax (%):', emax_var, 1, 200),
            ('EC50 (mg/L):', ec50_var, 0.01, 100),
            ('Hill n:', hill_var, 0.1, 10),
        ]
        for i, (label, var, lo, hi) in enumerate(params):
            row = i // 4
            col = i % 4
            tk.Label(params_grid, text=label,
                      bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                      font=('Arial', 9), anchor='e', width=12
                      ).grid(row=row, column=col*2, sticky='e', padx=2, pady=2)
            tk.Entry(params_grid, textvariable=var,
                      bg=Theme.c('panel2'), fg=Theme.c('text'),
                      insertbackground='#1c2128', relief='flat',
                      font=('Courier', 10), width=8,
                      highlightthickness=1, highlightbackground=Theme.c('border')
                      ).grid(row=row, column=col*2+1, sticky='w', padx=2, pady=2, ipady=2)

        # Route
        route_row = tk.Frame(ctrl_frame, bg=Theme.c('bg'))
        route_row.pack(fill='x', padx=14, pady=4)
        tk.Label(route_row, text=('Способ введения:' if self.lang == 'ru' else 'Route:'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(side='left')
        for val, lbl in [('oral', 'Oral'), ('iv', 'IV bolus')]:
            tk.Radiobutton(route_row, text=lbl, value=val,
                            variable=route_var,
                            bg=Theme.c('bg'), fg='#0550ae',
                            selectcolor='#eef1f5',
                            activebackground='#f5f7fa',
                            font=('Arial', 9)).pack(side='left', padx=8)

        # Results: plot + stats
        results_frame = tk.LabelFrame(dlg,
                                          text='  ' + ('Результаты PK/PD' if self.lang == 'ru'
                                                            else 'PK/PD Results') + '  ',
                                          bg=Theme.c('bg'), fg='#1a7f37',
                                          font=('Arial', 10, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=14, pady=4)

        plot_holder = tk.Frame(results_frame, bg=Theme.c('bg'))
        plot_holder.pack(fill='both', expand=True, padx=4, pady=4)
        stats_lbl = tk.Label(results_frame, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                                font=('Courier', 10), justify='left', anchor='w')
        stats_lbl.pack(fill='x', padx=14, pady=4, anchor='w')

        def run_pkpd():
            try:
                dose = dose_var.get()
                ka = ka_var.get()
                ke = ke_var.get()
                V = V_var.get()
                t_end = t_end_var.get()
                route = route_var.get()
                emax = emax_var.get()
                ec50 = ec50_var.get()
                hill = hill_var.get()
            except Exception as e:
                messagebox.showerror('VIREON', 'Bad parameter: ' + str(e))
                return

            # Time points
            n_pts = 100
            t_pts = [i * t_end / n_pts for i in range(n_pts + 1)]
            # PK
            conc = pk_one_compartment(dose, ka, ke, V, t_pts, route=route)
            # PD: Effect = Emax * C^n / (EC50^n + C^n)
            effect = [pd_emax_model(c, emax, ec50, hill) for c in conc]
            # Stats
            params_pk = pk_parameters(t_pts, conc)

            # Plot
            for w in plot_holder.winfo_children():
                w.destroy()
            try:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4),
                                                  facecolor='#f5f7fa')
                ax1.set_facecolor('#ffffff')
                ax2.set_facecolor('#ffffff')

                # PK
                ax1.plot(t_pts, conc, color='#8250df', linewidth=2,
                          label='Concentration')
                ax1.axhline(ec50, color='#9a6700', linestyle='--',
                              alpha=0.5, label='EC50')
                ax1.set_xlabel('Time (h)', color='#1c2128')
                ax1.set_ylabel('C (mg/L)', color='#1c2128')
                ax1.set_title('PK: Concentration vs Time', color='#0550ae')
                ax1.tick_params(colors='#57606a')
                ax1.grid(True, alpha=0.2, color='#d0d7de')
                ax1.legend(loc='upper right', facecolor='#ffffff',
                            edgecolor='#d0d7de', labelcolor='#1c2128')
                for sp in ax1.spines.values(): sp.set_color('#d0d7de')

                # PD
                ax2.plot(t_pts, effect, color='#1a7f37', linewidth=2,
                          label='Effect')
                ax2.set_xlabel('Time (h)', color='#1c2128')
                ax2.set_ylabel('Effect (%)', color='#1c2128')
                ax2.set_title('PD: Effect vs Time', color='#0550ae')
                ax2.tick_params(colors='#57606a')
                ax2.grid(True, alpha=0.2, color='#d0d7de')
                for sp in ax2.spines.values(): sp.set_color('#d0d7de')

                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=plot_holder)
                cv.get_tk_widget().pack(fill='both', expand=True)
            except Exception as e:
                biosim_log('warn', 'run_pkpd: {}'.format(e))
                tk.Label(plot_holder, text='Plot error: ' + str(e),
                          bg=Theme.c('bg'), fg='#cf222e').pack()

            # Stats display
            stats = ('PK Parameters:                              PD Parameters:\n'
                      '  Cmax    = {:.3f} mg/L                       Emax    = {:.1f}%\n'.format(
                          params_pk.get('Cmax', 0), emax) +
                      '  Tmax    = {:.2f} h                          EC50    = {:.2f} mg/L\n'.format(
                          params_pk.get('Tmax', 0), ec50) +
                      '  AUC     = {:.2f} mg·h/L                     Hill n  = {:.1f}\n'.format(
                          params_pk.get('AUC', 0), hill) +
                      '  t_1/2   = {:.2f} h                          Max effect: {:.1f}%\n'.format(
                          params_pk.get('t_half', 0), max(effect)))
            stats_lbl.config(text=stats)

        tk.Button(ctrl_frame, text='▶ ' + ('Запустить' if self.lang == 'ru' else 'Run'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 11, 'bold'), cursor='hand2',
                   padx=18, pady=4, command=run_pkpd
                   ).pack(pady=4)

        run_pkpd()  # initial

        tk.Button(dlg, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy).pack(fill='x', padx=14, pady=8)

    def _open_survival_analysis(self):
        """Kaplan-Meier survival analysis."""
        dlg = tk.Toplevel(self)
        dlg.title('Survival Analysis (Kaplan-Meier)')
        dlg.geometry('1000x720')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Анализ выживаемости' if self.lang == 'ru'
                                       else 'Survival Analysis'),
                  bg=Theme.c('bg'), fg='#cf222e',
                  font=('Arial', 15, 'bold')).pack(pady=10)

        # Two-group input
        groups_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        groups_frame.pack(fill='x', padx=14, pady=4)

        # Group A
        ga_frame = tk.LabelFrame(groups_frame,
                                     text='  Group A (Treatment)  ',
                                     bg=Theme.c('bg'), fg='#1a7f37',
                                     font=('Arial', 10, 'bold'), relief='flat')
        ga_frame.pack(side='left', fill='both', expand=True, padx=4)
        tk.Label(ga_frame,
                  text=('time,event (1=death,0=censored) per line:'
                        if self.lang == 'en' else
                        'время,событие (1=смерть,0=цензурировано) на строку:'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 9, 'italic')).pack(anchor='w', padx=4)
        ga_text = tk.Text(ga_frame, bg=Theme.c('panel2'), fg=Theme.c('text'),
                            insertbackground='#1c2128', relief='flat',
                            font=('Courier', 10), height=8, padx=8, pady=4)
        ga_text.pack(fill='both', expand=True, padx=4, pady=4)
        # Example data
        ga_text.insert('1.0',
            '5,1\n7,1\n10,0\n12,1\n15,1\n20,0\n22,1\n25,1\n30,0\n35,1')

        # Group B
        gb_frame = tk.LabelFrame(groups_frame,
                                     text='  Group B (Control)  ',
                                     bg=Theme.c('bg'), fg='#8250df',
                                     font=('Arial', 10, 'bold'), relief='flat')
        gb_frame.pack(side='left', fill='both', expand=True, padx=4)
        tk.Label(gb_frame,
                  text=('time,event per line:'
                        if self.lang == 'en' else 'время,событие на строку:'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 9, 'italic')).pack(anchor='w', padx=4)
        gb_text = tk.Text(gb_frame, bg=Theme.c('panel2'), fg=Theme.c('text'),
                            insertbackground='#1c2128', relief='flat',
                            font=('Courier', 10), height=8, padx=8, pady=4)
        gb_text.pack(fill='both', expand=True, padx=4, pady=4)
        gb_text.insert('1.0',
            '3,1\n5,1\n6,1\n8,1\n10,1\n11,1\n14,0\n16,1\n18,1\n20,1')

        # Results
        results_frame = tk.LabelFrame(dlg,
                                          text='  ' + ('K-M кривые' if self.lang == 'ru'
                                                            else 'K-M Curves') + '  ',
                                          bg=Theme.c('bg'), fg='#1a7f37',
                                          font=('Arial', 10, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=14, pady=4)

        plot_holder = tk.Frame(results_frame, bg=Theme.c('bg'))
        plot_holder.pack(fill='both', expand=True, padx=4, pady=4)
        stats_lbl = tk.Label(results_frame, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                                font=('Courier', 10), justify='left', anchor='w')
        stats_lbl.pack(fill='x', padx=14, pady=4, anchor='w')

        def parse_data(text):
            times, events = [], []
            for line in text.strip().split('\n'):
                line = line.strip()
                if not line or ',' not in line:
                    continue
                parts = line.split(',')
                try:
                    t = float(parts[0].strip())
                    e = int(parts[1].strip())
                    times.append(t)
                    events.append(e)
                except Exception as _exc:
                    biosim_log('warn', 'parse_data: {}'.format(_exc))
                    continue
            return times, events

        def run_km():
            ta, ea = parse_data(ga_text.get('1.0', 'end'))
            tb, eb = parse_data(gb_text.get('1.0', 'end'))
            if not ta or not tb:
                messagebox.showerror('VIREON', 'Need data in both groups')
                return

            km_a = kaplan_meier(ta, ea)
            km_b = kaplan_meier(tb, eb)
            lr = log_rank_test(ta, ea, tb, eb)

            # Plot
            for w in plot_holder.winfo_children():
                w.destroy()
            try:
                fig, ax = plt.subplots(figsize=(9, 4), facecolor='#f5f7fa')
                ax.set_facecolor('#ffffff')
                # Step plot (where='post')
                ax.step(km_a['time'], km_a['survival'], where='post',
                          color='#1a7f37', linewidth=2.5, label='Group A (n=' +
                          str(lr['n_a']) + ')')
                ax.step(km_b['time'], km_b['survival'], where='post',
                          color='#8250df', linewidth=2.5, label='Group B (n=' +
                          str(lr['n_b']) + ')')
                ax.set_xlabel('Time', color='#1c2128')
                ax.set_ylabel('Survival probability S(t)', color='#1c2128')
                ax.set_title('Kaplan-Meier Survival Curves', color='#0550ae')
                ax.tick_params(colors='#57606a')
                ax.grid(True, alpha=0.2, color='#d0d7de')
                ax.legend(loc='upper right', facecolor='#ffffff',
                            edgecolor='#d0d7de', labelcolor='#1c2128')
                ax.set_ylim([0, 1.05])
                for sp in ax.spines.values(): sp.set_color('#d0d7de')
                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=plot_holder)
                cv.get_tk_widget().pack(fill='both', expand=True)
            except Exception as e:
                biosim_log('warn', 'run_km: {}'.format(e))
                tk.Label(plot_holder, text='Plot error: ' + str(e),
                          bg=Theme.c('bg'), fg='#cf222e').pack()

            # Stats
            p = lr['p_value']
            sig = 'Significant' if p < 0.05 else '○ Not significant'
            color_sig = '#1a7f37' if p < 0.05 else '#9a6700'
            stats = ('Log-rank test:\n'
                      '  χ² = {:.3f}\n'.format(lr['chi_squared']) +
                      '  p = {:.4f}\n'.format(p) +
                      '  ' + sig + ' (α=0.05)\n')
            stats_lbl.config(text=stats, fg=color_sig)

        btn_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        btn_frame.pack(fill='x', padx=14, pady=4)
        tk.Button(btn_frame, text='▶ ' + ('Запустить K-M' if self.lang == 'ru'
                                              else 'Run K-M analysis'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 11, 'bold'), cursor='hand2',
                   padx=16, pady=4, command=run_km
                   ).pack(side='left', padx=2)
        tk.Button(btn_frame, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=14, pady=4, command=dlg.destroy
                   ).pack(side='right', padx=2)

        run_km()  # initial

    def _open_statistics(self):
        """Classical statistical tests dialog."""
        dlg = tk.Toplevel(self)
        dlg.title('Statistical Tests')
        dlg.geometry('900x700')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Статистические тесты' if self.lang == 'ru'
                                       else 'Statistical Tests'),
                  bg=Theme.c('bg'), fg='#0550ae',
                  font=('Arial', 15, 'bold')).pack(pady=10)

        # Test type selection
        test_var = tk.StringVar(value='t_test')
        test_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        test_frame.pack(fill='x', padx=14, pady=4)
        tk.Label(test_frame, text=('Тест:' if self.lang == 'ru' else 'Test:'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 10, 'bold')).pack(side='left', padx=4)
        for val, lbl in [
            ('t_test', "Student's t-test"),
            ('mann_whitney', 'Mann-Whitney U'),
            ('anova', 'One-way ANOVA'),
        ]:
            tk.Radiobutton(test_frame, text=lbl, value=val, variable=test_var,
                            bg=Theme.c('bg'), fg='#0550ae',
                            selectcolor='#eef1f5',
                            activebackground='#f5f7fa',
                            font=('Arial', 10, 'bold')
                            ).pack(side='left', padx=4)

        # Groups (3 groups, third only for ANOVA)
        groups_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        groups_frame.pack(fill='both', expand=True, padx=14, pady=4)

        group_texts = {}
        for i, (label, color) in enumerate([
            ('Group A', '#1a7f37'),
            ('Group B', '#8250df'),
            ('Group C (ANOVA only)', '#9a6700'),
        ]):
            gf = tk.LabelFrame(groups_frame, text='  ' + label + '  ',
                                 bg=Theme.c('bg'), fg=color,
                                 font=('Arial', 10, 'bold'), relief='flat')
            gf.pack(side='left', fill='both', expand=True, padx=4)
            tk.Label(gf, text=('Числа (по одному на строку):' if self.lang == 'ru'
                                  else 'Numbers (one per line):'),
                      bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                      font=('Arial', 9, 'italic')).pack(anchor='w', padx=4)
            txt = tk.Text(gf, bg=Theme.c('panel2'), fg=Theme.c('text'),
                            insertbackground='#1c2128', relief='flat',
                            font=('Courier', 10), height=10, padx=8, pady=4)
            txt.pack(fill='both', expand=True, padx=4, pady=4)
            group_texts[i] = txt
        # Sample data
        group_texts[0].insert('1.0', '10\n12\n11\n13\n15\n10\n11\n14\n12\n13')
        group_texts[1].insert('1.0', '15\n17\n16\n18\n20\n14\n16\n19\n18\n17')
        group_texts[2].insert('1.0', '20\n22\n21\n23\n25\n19\n21\n24\n22\n23')

        # Results
        results_lbl = tk.Label(dlg, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                                  font=('Courier', 11), justify='left', anchor='w')
        results_lbl.pack(fill='x', padx=20, pady=8, anchor='w')

        def parse_vals(text):
            vals = []
            for line in text.strip().split('\n'):
                try:
                    vals.append(float(line.strip()))
                except Exception as _exc:
                    biosim_log('warn', 'parse_vals: {}'.format(_exc))
                    continue
            return vals

        def run_test():
            test = test_var.get()
            ga = parse_vals(group_texts[0].get('1.0', 'end'))
            gb = parse_vals(group_texts[1].get('1.0', 'end'))
            gc = parse_vals(group_texts[2].get('1.0', 'end'))

            if test == 't_test':
                if len(ga) < 2 or len(gb) < 2:
                    results_lbl.config(text='Need ≥2 values per group', fg='#cf222e')
                    return
                r = t_test(ga, gb)
                if not r:
                    return
                p = r['p_value']
                sig = 'Significant' if p < 0.05 else '○ Not significant'
                color = '#1a7f37' if p < 0.05 else '#9a6700'
                ci_a = confidence_interval_95(ga)
                ci_b = confidence_interval_95(gb)
                txt = ("Student's t-test (Welch):\n"
                        '  t = {:.4f}, df = {:.2f}\n'.format(r['t'], r['df']) +
                        '  p = {:.5f}\n\n'.format(p) +
                        '  Group A: mean = {:.3f} ± {:.3f} (n={})\n'.format(
                            r['mean_a'], r['std_a'], r['n_a']) +
                        '           95% CI: [{:.3f}, {:.3f}]\n'.format(
                            ci_a['lower'], ci_a['upper']) +
                        '  Group B: mean = {:.3f} ± {:.3f} (n={})\n'.format(
                            r['mean_b'], r['std_b'], r['n_b']) +
                        '           95% CI: [{:.3f}, {:.3f}]\n\n'.format(
                            ci_b['lower'], ci_b['upper']) +
                        '  Difference: {:.3f} (SE = {:.3f})\n\n'.format(
                            r['diff'], r['se']) +
                        '  ' + sig + ' (α = 0.05)')
                results_lbl.config(text=txt, fg=color)

            elif test == 'mann_whitney':
                if len(ga) < 2 or len(gb) < 2:
                    results_lbl.config(text='Need ≥2 values per group', fg='#cf222e')
                    return
                r = mann_whitney_u(ga, gb)
                if not r:
                    return
                p = r['p_value']
                sig = 'Significant' if p < 0.05 else '○ Not significant'
                color = '#1a7f37' if p < 0.05 else '#9a6700'
                txt = ('Mann-Whitney U test:\n'
                        '  U = {:.2f}\n'.format(r['U']) +
                        '  z = {:.4f}\n'.format(r['z']) +
                        '  p = {:.5f}\n\n'.format(p) +
                        '  n_a = {}, n_b = {}\n\n'.format(r['n_a'], r['n_b']) +
                        '  ' + sig + ' (α = 0.05)\n\n'
                        '  Non-parametric — does not assume normality.')
                results_lbl.config(text=txt, fg=color)

            elif test == 'anova':
                if len(ga) < 2 or len(gb) < 2:
                    results_lbl.config(text='Need ≥2 values per group', fg='#cf222e')
                    return
                # Use all 3 groups if C non-empty, else 2
                if len(gc) >= 2:
                    r = anova_oneway(ga, gb, gc)
                else:
                    r = anova_oneway(ga, gb)
                if not r:
                    return
                p = r['p_value']
                sig = 'Significant' if p < 0.05 else '○ Not significant'
                color = '#1a7f37' if p < 0.05 else '#9a6700'
                txt = ('One-way ANOVA (F-test):\n'
                        '  F = {:.4f}\n'.format(r['F']) +
                        '  df = ({}, {})\n'.format(r['df_between'], r['df_within']) +
                        '  p ≈ {:.5f}\n\n'.format(p) +
                        '  Groups: k = {}, total N = {}\n'.format(r['k'], r['N']) +
                        '  SS between: {:.3f}\n'.format(r['SSB']) +
                        '  SS within:  {:.3f}\n\n'.format(r['SSW']) +
                        '  ' + sig + ' (α = 0.05)\n\n'
                        '  Note: post-hoc test needed to find which groups differ.')
                results_lbl.config(text=txt, fg=color)

        # Bottom
        bottom = tk.Frame(dlg, bg=Theme.c('bg'))
        bottom.pack(fill='x', padx=14, pady=8)
        tk.Button(bottom, text='▶ ' + ('Запустить' if self.lang == 'ru' else 'Run test'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 11, 'bold'), cursor='hand2',
                   padx=18, pady=4, command=run_test
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=14, pady=4, command=dlg.destroy
                   ).pack(side='right', padx=2)

        run_test()  # initial

    def _open_data_import(self):
        """Import patient cohort data from CSV/Excel for analysis."""
        from tkinter import filedialog
        filepath = filedialog.askopenfilename(
            title='Import data',
            filetypes=[('CSV', '*.csv'), ('Excel', '*.xlsx *.xls'),
                        ('TSV', '*.tsv'), ('All', '*.*')])
        if not filepath:
            return

        try:
            import pandas as pd
            if filepath.endswith(('.xlsx', '.xls')):
                df = pd.read_excel(filepath)
            elif filepath.endswith('.tsv'):
                df = pd.read_csv(filepath, sep='\t')
            else:
                df = pd.read_csv(filepath)
        except ImportError:
            messagebox.showerror('VIREON',
                'pandas not installed. Install: pip install pandas openpyxl')
            return
        except Exception as e:
            messagebox.showerror('VIREON', 'Read error: ' + str(e))
            return

        # Show data
        dlg = tk.Toplevel(self)
        dlg.title('Imported Data')
        dlg.geometry('960x600')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text='' + os.path.basename(filepath),
                  bg=Theme.c('bg'), fg='#1a7f37',
                  font=('Arial', 14, 'bold')).pack(pady=10)
        tk.Label(dlg,
                  text='{} rows × {} columns'.format(len(df), len(df.columns)),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 10, 'italic')).pack()

        # Show as table (tkinter Treeview)
        from tkinter import ttk as _ttk
        tree_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        tree_frame.pack(fill='both', expand=True, padx=14, pady=8)

        sb_v = tk.Scrollbar(tree_frame, orient='vertical')
        sb_v.pack(side='right', fill='y')
        sb_h = tk.Scrollbar(tree_frame, orient='horizontal')
        sb_h.pack(side='bottom', fill='x')

        style = _ttk.Style()
        try:
            style.configure('Imported.Treeview',
                              background='#ffffff', foreground=Theme.c('text'),
                              fieldbackground='#ffffff',
                              borderwidth=0)
            style.configure('Imported.Treeview.Heading',
                              background='#eef1f5', foreground='#1a7f37')
        except Exception as _exc:
            biosim_log('warn', '_open_data_import: {}'.format(_exc))

        tree = _ttk.Treeview(tree_frame, style='Imported.Treeview',
                                yscrollcommand=sb_v.set,
                                xscrollcommand=sb_h.set)
        tree.pack(fill='both', expand=True)
        sb_v.config(command=tree.yview)
        sb_h.config(command=tree.xview)

        columns = list(df.columns)
        tree['columns'] = columns
        tree['show'] = 'headings'
        for col in columns:
            tree.heading(col, text=col)
            tree.column(col, width=120, anchor='w')

        # Insert rows (limit to 1000)
        for idx, row in df.head(1000).iterrows():
            tree.insert('', 'end', values=list(row))

        # Summary stats
        info_text = ('Columns: ' + ', '.join(columns[:5]) +
                      ('...' if len(columns) > 5 else '') + '\n')
        try:
            numeric_cols = df.select_dtypes(include='number').columns
            for c in numeric_cols[:5]:
                info_text += '  {}: mean={:.2f}, std={:.2f}\n'.format(
                    c, df[c].mean(), df[c].std())
        except Exception as _exc:
            biosim_log('warn', '_open_data_import: {}'.format(_exc))
        tk.Label(dlg, text=info_text, bg=Theme.c('bg'), fg=Theme.c('text'),
                  font=('Courier', 9), justify='left', anchor='w'
                  ).pack(padx=14, pady=4, anchor='w')

        # Save to clipboard for stats dialog
        self._imported_data = df

        bottom = tk.Frame(dlg, bg=Theme.c('bg'))
        bottom.pack(fill='x', padx=14, pady=8)
        tk.Label(bottom,
                  text=('Данные доступны для статистики и графиков.' if self.lang == 'ru'
                        else 'Data available for stats and plots.'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 9, 'italic')).pack(side='left')
        tk.Button(bottom, text='Close', bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                   relief='solid', bd=1, font=('Arial', 10), cursor='hand2',
                   padx=14, pady=4, command=dlg.destroy
                   ).pack(side='right', padx=2)

    def _open_report_generator(self):
        """Generate scientific report based on current state."""
        if not self.cur_dis:
            messagebox.showinfo('VIREON',
                ('Сначала выбери болезнь и запусти симуляцию!' if self.lang == 'ru'
                  else 'Choose a disease and run simulation first!'))
            return

        dlg = tk.Toplevel(self)
        dlg.title('Scientific Report Generator')
        dlg.geometry('900x680')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Генератор научного отчёта'
                                       if self.lang == 'ru' else
                                       'Scientific Report Generator'),
                  bg=Theme.c('bg'), fg='#8250df',
                  font=('Arial', 15, 'bold')).pack(pady=10)

        # Format selection
        fmt_var = tk.StringVar(value='markdown')
        fmt_row = tk.Frame(dlg, bg=Theme.c('bg'))
        fmt_row.pack(fill='x', padx=14, pady=4)
        tk.Label(fmt_row, text=('Формат:' if self.lang == 'ru' else 'Format:'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 10, 'bold')).pack(side='left', padx=4)
        for val, lbl in [('markdown', 'Markdown'),
                          ('latex', 'LaTeX'),
                          ('html', 'HTML')]:
            tk.Radiobutton(fmt_row, text=lbl, value=val, variable=fmt_var,
                            bg=Theme.c('bg'), fg='#0550ae',
                            selectcolor='#eef1f5',
                            activebackground='#f5f7fa',
                            font=('Arial', 10, 'bold')
                            ).pack(side='left', padx=4)

        # Preview
        preview_frame = tk.LabelFrame(dlg,
                                          text='  ' + ('Превью' if self.lang == 'ru'
                                                            else 'Preview') + '  ',
                                          bg=Theme.c('bg'), fg='#1a7f37',
                                          font=('Arial', 10, 'bold'), relief='flat')
        preview_frame.pack(fill='both', expand=True, padx=14, pady=4)

        sb = tk.Scrollbar(preview_frame)
        sb.pack(side='right', fill='y')
        preview_txt = tk.Text(preview_frame, bg=Theme.c('panel'), fg=Theme.c('text'),
                                  font=('Courier', 10), relief='flat',
                                  wrap='word', padx=10, pady=8,
                                  yscrollcommand=sb.set)
        preview_txt.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        sb.config(command=preview_txt.yview)

        def generate_report():
            fmt = fmt_var.get()
            content = self._build_scientific_report(fmt)
            preview_txt.delete('1.0', 'end')
            preview_txt.insert('1.0', content)

        def save_report():
            fmt = fmt_var.get()
            from tkinter import filedialog
            ext_map = {'markdown': '.md', 'latex': '.tex', 'html': '.html'}
            ext = ext_map.get(fmt, '.txt')
            filepath = filedialog.asksaveasfilename(
                title='Save report',
                defaultextension=ext,
                initialfile='VIREON_Report_' + (self.cur_dis or 'report') + ext,
                filetypes=[(fmt.upper(), '*' + ext), ('All', '*.*')])
            if not filepath:
                return
            try:
                with open(filepath, 'w', encoding='utf-8') as f:
                    f.write(preview_txt.get('1.0', 'end-1c'))
                messagebox.showinfo('VIREON', 'Saved: ' + filepath)
            except Exception as e:
                messagebox.showerror('VIREON', str(e))

        # Generate on format change
        fmt_var.trace('w', lambda *a: generate_report())
        generate_report()

        # Bottom
        bottom = tk.Frame(dlg, bg=Theme.c('bg'))
        bottom.pack(fill='x', padx=14, pady=8)
        tk.Button(bottom, text='Regenerate',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=12, pady=4, command=generate_report
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text='Save',
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=save_report
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=14, pady=4, command=dlg.destroy
                   ).pack(side='right', padx=2)

    def _build_scientific_report(self, fmt):
        """Generate report content in requested format."""
        import datetime
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
        disease = self.cur_dis or 'N/A'
        try:
            disease_name = (self.cur_dat.get('name_' + self.lang)
                              or self.cur_dat.get('name_ru', disease))
        except Exception as _exc:
            biosim_log('warn', '_build_scientific_report: {}'.format(_exc))
            disease_name = disease
        drug = self.drug_var.get() or 'N/A'
        try:
            dose = self.dose_mg_var.get()
        except Exception as _exc:
            biosim_log('warn', '_build_scientific_report: {}'.format(_exc))
            dose = 'N/A'
        try:
            efficacy = self.drug_dose.get()
        except Exception as _exc:
            biosim_log('warn', '_build_scientific_report: {}'.format(_exc))
            efficacy = 'N/A'

        # Profile
        try:
            prof = self._load_patient_profile()
        except Exception as _exc:
            biosim_log('warn', '_build_scientific_report: {}'.format(_exc))
            prof = {}

        if fmt == 'markdown':
            out = ['# VIREON Scientific Report\n',
                    '*Generated: ' + now + ' · VIREON v5.0*\n',
                    '\n## 1. Subject Information\n',
                    '| Attribute | Value |',
                    '|-----------|-------|',
                    '| Patient ID | ' + (prof.get('name') or 'Anonymous') + ' |',
                    '| Age | ' + str(prof.get('age', 'N/A')) + ' |',
                    '| Sex | ' + str(prof.get('sex', 'N/A')) + ' |',
                    '| Weight | ' + str(prof.get('weight', 'N/A')) + ' kg |',
                    '| BMI | ' + str(prof.get('bmi', 'N/A')) + ' |',
                    '\n## 2. Disease\n',
                    '**Diagnosis:** ' + disease_name + ' (' + disease + ')\n',
                    '**Affected organs:** ' + str(self.cur_dat.get('organs', 'N/A') if self.cur_dat else 'N/A') + '\n',
                    '**Symptoms:** ' + str(self.cur_dat.get('symptoms_' + self.lang, 'N/A')
                                              if self.cur_dat else 'N/A') + '\n',
                    '\n## 3. Treatment Protocol\n',
                    '| Parameter | Value |',
                    '|-----------|-------|',
                    '| Drug | ' + drug + ' |',
                    '| Dose | ' + str(dose) + ' mg |',
                    '| Efficacy parameter | ' + str(efficacy) + '% |',
                    '\n## 4. Methods\n',
                    'Simulation performed using VIREON v5.0 with built-in disease '
                    'progression models. Drug effect modeled via dose-response curve. '
                    'For details see https://github.com/biosim\n',
                    '\n## 5. Results\n',
                    '*Simulation results to be added after running.*\n',
                    '\n## 6. Discussion\n',
                    '*To be filled by researcher.*\n',
                    '\n## 7. References\n',
                    '1. VIREON Software (v5.0). https://github.com/biosim\n',
                    '\n---\n',
                    '*Report generated automatically. For research/educational use only. '
                    'Not a substitute for medical advice.*']
            return '\n'.join(out)

        elif fmt == 'latex':
            out = [r'\documentclass{article}',
                    r'\usepackage[utf8]{inputenc}',
                    r'\usepackage{booktabs}',
                    r'\usepackage{hyperref}',
                    r'\title{VIREON Scientific Report}',
                    r'\author{VIREON v5.0}',
                    r'\date{' + now + '}',
                    r'\begin{document}',
                    r'\maketitle',
                    r'\section{Subject Information}',
                    r'\begin{tabular}{ll}',
                    r'\toprule',
                    r'Attribute & Value \\',
                    r'\midrule',
                    'Patient ID & ' + (prof.get('name') or 'Anonymous') + r' \\',
                    'Age & ' + str(prof.get('age', 'N/A')) + r' \\',
                    'Sex & ' + str(prof.get('sex', 'N/A')) + r' \\',
                    'Weight & ' + str(prof.get('weight', 'N/A')) + r' kg \\',
                    'BMI & ' + str(prof.get('bmi', 'N/A')) + r' \\',
                    r'\bottomrule',
                    r'\end{tabular}',
                    r'\section{Disease}',
                    r'\textbf{Diagnosis:} ' + disease_name + ' (' + disease + r')',
                    r'\section{Treatment Protocol}',
                    r'\begin{itemize}',
                    r'  \item Drug: ' + drug,
                    r'  \item Dose: ' + str(dose) + ' mg',
                    r'  \item Efficacy: ' + str(efficacy) + r'\%',
                    r'\end{itemize}',
                    r'\section{Methods}',
                    'Simulation performed using VIREON v5.0.',
                    r'\section{Results}',
                    r'\textit{To be added.}',
                    r'\section{Discussion}',
                    r'\textit{To be filled by researcher.}',
                    r'\section{References}',
                    r'\begin{enumerate}',
                    r'  \item VIREON Software (v5.0). \url{https://github.com/biosim}',
                    r'\end{enumerate}',
                    r'\end{document}']
            return '\n'.join(out)

        else:  # html
            out = ['<!DOCTYPE html><html><head><meta charset="UTF-8">',
                    '<title>VIREON Report</title>',
                    '<style>body{font-family:Arial,sans-serif;max-width:800px;'
                    'margin:0 auto;padding:30px;color:#222;}'
                    'h1{color:#0969da;border-bottom:2px solid #ddd;}'
                    'h2{color:#1a7f37;margin-top:24px;}'
                    'table{border-collapse:collapse;width:100%;}'
                    'th,td{border:1px solid #ddd;padding:8px;text-align:left;}'
                    'th{background:#f0f0f0;}</style></head><body>',
                    '<h1>VIREON Scientific Report</h1>',
                    '<p><em>Generated: ' + now + ' · VIREON v5.0</em></p>',
                    '<h2>1. Subject Information</h2>',
                    '<table>',
                    '<tr><th>Attribute</th><th>Value</th></tr>',
                    '<tr><td>Patient ID</td><td>' + (prof.get('name') or 'Anonymous') + '</td></tr>',
                    '<tr><td>Age</td><td>' + str(prof.get('age', 'N/A')) + '</td></tr>',
                    '<tr><td>Sex</td><td>' + str(prof.get('sex', 'N/A')) + '</td></tr>',
                    '<tr><td>Weight</td><td>' + str(prof.get('weight', 'N/A')) + ' kg</td></tr>',
                    '<tr><td>BMI</td><td>' + str(prof.get('bmi', 'N/A')) + '</td></tr>',
                    '</table>',
                    '<h2>2. Disease</h2>',
                    '<p><strong>Diagnosis:</strong> ' + disease_name + '</p>',
                    '<h2>3. Treatment Protocol</h2>',
                    '<ul>',
                    '<li>Drug: ' + drug + '</li>',
                    '<li>Dose: ' + str(dose) + ' mg</li>',
                    '<li>Efficacy: ' + str(efficacy) + '%</li>',
                    '</ul>',
                    '<h2>4. Methods</h2>',
                    '<p>Simulation performed using VIREON v5.0.</p>',
                    '<h2>5. Results</h2>',
                    '<p><em>To be added.</em></p>',
                    '<h2>6. Discussion</h2>',
                    '<p><em>To be filled by researcher.</em></p>',
                    '<h2>7. References</h2>',
                    '<ol><li>VIREON Software (v5.0)</li></ol>',
                    '<hr><p style="color:#888;font-size:0.9em">'
                    'For research/educational use only.</p>',
                    '</body></html>']
            return '\n'.join(out)

    def _open_manhattan_plot(self):
        """Manhattan plot for GWAS-style data."""
        dlg = tk.Toplevel(self)
        dlg.title('Manhattan Plot (GWAS-style)')
        dlg.geometry('1100x700')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        tk.Label(dlg, text='Manhattan Plot — ' + ('GWAS визуализация'
                                                         if self.lang == 'ru' else
                                                         'GWAS Visualization'),
                  bg=Theme.c('bg'), fg='#1a7f37',
                  font=('Arial', 14, 'bold')).pack(pady=10)
        tk.Label(dlg,
                  text=('Загрузи CSV: columns = chromosome, position, p_value'
                        if self.lang == 'en' else
                        'Загрузи CSV: колонки = chromosome, position, p_value'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                  font=('Arial', 10, 'italic')).pack()

        # Plot area
        plot_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        plot_frame.pack(fill='both', expand=True, padx=14, pady=8)
        plot_holder = tk.Frame(plot_frame, bg=Theme.c('panel'), height=500)
        plot_holder.pack(fill='both', expand=True, padx=4, pady=4)

        # Stats
        stats_lbl = tk.Label(plot_frame, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                                font=('Courier', 10), justify='left')
        stats_lbl.pack(fill='x', padx=14, pady=2)

        def draw_manhattan(df=None):
            for w in plot_holder.winfo_children():
                w.destroy()
            try:
                if df is None:
                    # Generate synthetic example
                    import random
                    import math as _m
                    rng = random.Random(42)
                    chroms = []
                    poss = []
                    pvals = []
                    for chrom in range(1, 23):
                        n_snps = rng.randint(50, 200)
                        for _ in range(n_snps):
                            chroms.append(chrom)
                            poss.append(rng.uniform(0, 250_000_000))
                            # Most SNPs not significant
                            if rng.random() < 0.005:
                                pvals.append(rng.uniform(1e-12, 1e-8))  # significant
                            elif rng.random() < 0.05:
                                pvals.append(rng.uniform(1e-7, 1e-5))  # suggestive
                            else:
                                pvals.append(rng.uniform(0.01, 1.0))   # not significant
                    df_data = list(zip(chroms, poss, pvals))
                else:
                    cols = list(df.columns)
                    chrom_col = next((c for c in cols if c.lower() in
                                         ['chromosome', 'chrom', 'chr']), cols[0])
                    pos_col = next((c for c in cols if c.lower() in
                                       ['position', 'pos', 'bp']), cols[1])
                    p_col = next((c for c in cols if 'p' in c.lower()
                                     and 'value' in c.lower()),
                                    next((c for c in cols if 'p' in c.lower()),
                                          cols[2]))
                    df_data = [(row[chrom_col], row[pos_col], row[p_col])
                                 for _, row in df.iterrows()]

                fig, ax = plt.subplots(figsize=(11, 5), facecolor='#f5f7fa')
                ax.set_facecolor('#ffffff')

                # Group by chromosome, compute cumulative position
                df_data.sort(key=lambda x: (x[0], x[1]))
                cum_pos = 0
                tick_pos = []
                tick_labels = []
                xs, ys, colors = [], [], []
                colors_alt = ['#0969da', '#8250df']
                current_chr = None
                chr_start = 0
                import math as _m
                for chrom, pos, pval in df_data:
                    if chrom != current_chr:
                        if current_chr is not None:
                            tick_pos.append((chr_start + cum_pos) / 2)
                            tick_labels.append(str(current_chr))
                        current_chr = chrom
                        chr_start = cum_pos
                    xs.append(cum_pos + pos / 1e6)
                    ys.append(-_m.log10(max(pval, 1e-300)))
                    colors.append(colors_alt[int(chrom) % 2 if isinstance(chrom, (int, float))
                                                  else hash(str(chrom)) % 2])
                    cum_pos = max(cum_pos, xs[-1])
                if current_chr is not None:
                    tick_pos.append((chr_start + cum_pos) / 2)
                    tick_labels.append(str(current_chr))

                ax.scatter(xs, ys, c=colors, s=3, alpha=0.7)
                # Significance lines
                ax.axhline(-_m.log10(5e-8), color='#cf222e',
                              linestyle='--', linewidth=1, label='Genome-wide (5×10⁻⁸)')
                ax.axhline(-_m.log10(1e-5), color='#9a6700',
                              linestyle=':', linewidth=1, label='Suggestive (1×10⁻⁵)')
                ax.set_xlabel('Chromosome', color='#1c2128')
                ax.set_ylabel('-log₁₀(p-value)', color='#1c2128')
                ax.set_title('Manhattan Plot', color='#0550ae')
                ax.set_xticks(tick_pos)
                ax.set_xticklabels(tick_labels, rotation=0, fontsize=8)
                ax.tick_params(colors='#57606a')
                ax.grid(True, alpha=0.2, color='#d0d7de', axis='y')
                ax.legend(loc='upper right', facecolor='#ffffff',
                            edgecolor='#d0d7de', labelcolor='#1c2128')
                for sp in ax.spines.values(): sp.set_color('#d0d7de')
                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=plot_holder)
                cv.get_tk_widget().pack(fill='both', expand=True)

                # Count significant SNPs
                n_sig = sum(1 for _, _, p in df_data if p < 5e-8)
                n_sug = sum(1 for _, _, p in df_data if 5e-8 <= p < 1e-5)
                stats_lbl.config(text='Total SNPs: ' + str(len(df_data)) +
                                       ' · Genome-wide significant (p<5×10⁻⁸): ' +
                                       str(n_sig) +
                                       ' · Suggestive (p<1×10⁻⁵): ' + str(n_sug))
            except Exception as e:
                biosim_log('warn', 'draw_manhattan: {}'.format(e))
                tk.Label(plot_holder, text='Plot error: ' + str(e),
                          bg=Theme.c('panel'), fg='#cf222e').pack(padx=20, pady=20)

        def load_csv():
            from tkinter import filedialog
            filepath = filedialog.askopenfilename(
                title='Load GWAS data',
                filetypes=[('CSV', '*.csv'), ('TSV', '*.tsv'), ('All', '*.*')])
            if not filepath:
                return
            try:
                import pandas as pd
                if filepath.endswith('.tsv'):
                    df = pd.read_csv(filepath, sep='\t')
                else:
                    df = pd.read_csv(filepath)
                draw_manhattan(df)
            except Exception as e:
                messagebox.showerror('VIREON', str(e))

        # Bottom
        bottom = tk.Frame(dlg, bg=Theme.c('bg'))
        bottom.pack(fill='x', padx=14, pady=8)
        tk.Button(bottom, text='Load CSV',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=load_csv
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text='Example data',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=lambda: draw_manhattan(None)
                   ).pack(side='left', padx=2)
        tk.Button(bottom, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=14, pady=4, command=dlg.destroy
                   ).pack(side='right', padx=2)

        draw_manhattan(None)  # show example

    def _autosave_tick(self):
        """Auto-save every 5 minutes."""
        try:
            import json
            settings_path = os.path.expanduser('~/.biosim/settings.json')
            do_autosave = True
            if os.path.exists(settings_path):
                with open(settings_path) as f:
                    s = json.load(f)
                do_autosave = s.get('autosave', True)
            if do_autosave and self.cur_dis:
                # Save autosave file
                path = os.path.expanduser('~/.biosim/autosave.biosim')
                os.makedirs(os.path.dirname(path), exist_ok=True)
                state = {
                    'disease': self.cur_dis,
                    'drug': self.drug_var.get(),
                    'dose_mg': self.dose_mg_var.get(),
                    'efficacy': self.drug_dose.get(),
                    'timestamp': str(__import__('datetime').datetime.now()),
                }
                with open(path, 'w') as f:
                    json.dump(state, f, indent=2, ensure_ascii=False)
        except Exception as e:
            biosim_log('error', 'autosave: ' + str(e))
        # Schedule next
        self.after(5 * 60 * 1000, self._autosave_tick)


    def _open_epidemiology(self):
        """Epidemiology simulator — SIR/SEIR models with interventions."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Эпидемиология' if is_ru else 'Epidemiology'))
        dlg.geometry('900x650')
        dlg.configure(bg=Theme.c('bg'))

        tk.Label(dlg, text=('Моделирование эпидемий (SIR/SEIR)'
                                       if is_ru else 'Epidemic modeling (SIR/SEIR)'),
                  bg=Theme.c('bg'), fg='#cf222e',
                  font=('Arial', 15, 'bold')).pack(pady=10)

        ctrl = tk.Frame(dlg, bg=Theme.c('bg')); ctrl.pack(fill='x', padx=16)

        # Parameters
        params = {}
        param_defs = [
            ('pop', 'Население' if is_ru else 'Population', 10000),
            ('I0', 'Заражённых вначале' if is_ru else 'Initial infected', 10),
            ('beta', 'β (заразность)' if is_ru else 'β (transmission)', 0.4),
            ('gamma', 'γ (выздоровление)' if is_ru else 'γ (recovery)', 0.1),
            ('days', 'Дней' if is_ru else 'Days', 120),
            ('vacc', 'Вакцинация %' if is_ru else 'Vaccination %', 0),
        ]
        for i, (key, label, default) in enumerate(param_defs):
            row = tk.Frame(ctrl, bg=Theme.c('bg'))
            row.grid(row=i // 3, column=i % 3, sticky='w', padx=8, pady=4)
            tk.Label(row, text=label + ':', bg=Theme.c('bg'), fg=Theme.c('text'),
                      font=('Arial', 9), width=18, anchor='w').pack(side='left')
            var = tk.StringVar(value=str(default))
            tk.Entry(row, textvariable=var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                      width=8, relief='flat').pack(side='left')
            params[key] = var

        model_var = tk.StringVar(value='SIR')
        mrow = tk.Frame(ctrl, bg=Theme.c('bg')); mrow.grid(row=2, column=0, columnspan=3, pady=6)
        tk.Label(mrow, text=('Модель:' if is_ru else 'Model:'), bg=Theme.c('bg'),
                  fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
        for m in ['SIR', 'SEIR']:
            tk.Radiobutton(mrow, text=m, variable=model_var, value=m,
                            bg=Theme.c('bg'), fg='#0969da', selectcolor='#eef1f5',
                            font=('Arial', 9)).pack(side='left', padx=6)

        # Result area with matplotlib
        fig_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        fig_frame.pack(fill='both', expand=True, padx=16, pady=8)

        info_lbl = tk.Label(dlg, text='', bg=Theme.c('bg'), fg='#1a7f37',
                             font=('Arial', 10, 'bold'))
        info_lbl.pack(pady=4)

        canvas_holder = {'canvas': None}

        def run_sim():
            try:
                pop = int(float(params['pop'].get()))
                I0 = int(float(params['I0'].get()))
                beta = float(params['beta'].get())
                gamma = float(params['gamma'].get())
                days = int(float(params['days'].get()))
                vacc = float(params['vacc'].get()) / 100.0
            except ValueError as _exc:
                biosim_log('warn', 'run_sim: {}'.format(_exc))
                info_lbl.config(text=('Проверьте числа' if is_ru else 'Check numbers'),
                                 fg='#cf222e')
                return

            if model_var.get() == 'SEIR':
                r = seir_model(pop - I0, I0, I0, 0, beta, 0.2, gamma, days)
                has_e = True
            else:
                if vacc > 0:
                    r = sir_with_vaccination(pop, I0, beta, gamma, vacc, days)
                else:
                    r = sir_model(pop - I0, I0, 0, beta, gamma, days)
                has_e = False

            # Plot
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if canvas_holder['canvas']:
                canvas_holder['canvas'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['S'], color='#0969da', label='S ' + ('восприимч.' if is_ru else 'Susceptible'))
            if has_e:
                ax.plot(r['times'], r['E'], color='#9a6700', label='E ' + ('латентн.' if is_ru else 'Exposed'))
            ax.plot(r['times'], r['I'], color='#cf222e', label='I ' + ('заражённ.' if is_ru else 'Infected'))
            ax.plot(r['times'], r['R'], color='#1a7f37', label='R ' + ('выздоров.' if is_ru else 'Recovered'))
            ax.set_xlabel('Дни' if is_ru else 'Days', color='#1c2128')
            ax.set_ylabel('Люди' if is_ru else 'People', color='#1c2128')
            ax.tick_params(colors='#57606a')
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            ax.grid(alpha=0.2)
            for spine in ax.spines.values():
                spine.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=fig_frame)
            cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True)
            canvas_holder['canvas'] = cv

            # Info
            txt = 'R₀ = {} | '.format(r['R0'])
            txt += ('Пик: ' if is_ru else 'Peak: ') + '{} '.format(r['peak_infected'])
            txt += ('на день ' if is_ru else 'on day ') + '{} | '.format(r['peak_day'])
            txt += ('Охват: ' if is_ru else 'Attack rate: ') + '{}%'.format(r['attack_rate'])
            if 'herd_immunity_threshold' in r:
                txt += ' | ' + ('Порог иммунитета: ' if is_ru else 'Herd threshold: ')
                txt += '{}%'.format(r['herd_immunity_threshold'])
            info_lbl.config(text=txt, fg='#1a7f37')

        tk.Button(ctrl, text=('▶ Запустить' if is_ru else '▶ Run'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run_sim).grid(row=3, column=0, columnspan=3,
                                          sticky='ew', pady=8)
        run_sim()

    def _open_neuroscience(self):
        """Neuroscience — Hodgkin-Huxley action potential simulator."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Нейробиология' if is_ru else 'Neuroscience'))
        dlg.geometry('900x650')
        dlg.configure(bg=Theme.c('bg'))

        tk.Label(dlg, text=('Нейрон Ходжкина-Хаксли'
                                       if is_ru else 'Hodgkin-Huxley neuron'),
                  bg=Theme.c('bg'), fg='#8250df',
                  font=('Arial', 15, 'bold')).pack(pady=10)
        tk.Label(dlg, text=('Модель потенциала действия нервной клетки'
                              if is_ru else 'Action potential model of a nerve cell'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack()

        ctrl = tk.Frame(dlg, bg=Theme.c('bg')); ctrl.pack(fill='x', padx=16, pady=8)
        tk.Label(ctrl, text=('Ток I (μA/cm²):' if is_ru else 'Current I (μA/cm²):'),
                  bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
        cur_var = tk.DoubleVar(value=10)
        cur_lbl = tk.Label(ctrl, text='10', bg=Theme.c('bg'), fg='#8250df',
                            font=('Arial', 10, 'bold'), width=5)
        cur_lbl.pack(side='right')
        tk.Scale(ctrl, from_=0, to=50, orient='horizontal', variable=cur_var,
                  bg=Theme.c('bg'), fg=Theme.c('text'), highlightthickness=0,
                  troughcolor='#eef1f5', resolution=1, showvalue=False,
                  command=lambda v: cur_lbl.config(text=str(int(float(v))))
                  ).pack(side='left', fill='x', expand=True, padx=8)

        fig_frame = tk.Frame(dlg, bg=Theme.c('bg'))
        fig_frame.pack(fill='both', expand=True, padx=16, pady=8)
        info_lbl = tk.Label(dlg, text='', bg=Theme.c('bg'), fg='#1a7f37',
                             font=('Arial', 10, 'bold'))
        info_lbl.pack(pady=4)
        ch = {'canvas': None}

        def run():
            I = cur_var.get()
            r = hodgkin_huxley(I_inj=I, t_max=50)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch['canvas']:
                ch['canvas'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4.2), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['V'], color='#8250df', linewidth=1.2)
            ax.set_xlabel(('Время (мс)' if is_ru else 'Time (ms)'), color='#1c2128')
            ax.set_ylabel(('Потенциал V (мВ)' if is_ru else 'Potential V (mV)'), color='#1c2128')
            ax.axhline(y=-65, color='#57606a', linestyle='--', alpha=0.4)
            ax.tick_params(colors='#57606a')
            ax.grid(alpha=0.2)
            for sp in ax.spines.values():
                sp.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=fig_frame)
            cv.draw(); cv.get_tk_widget().pack(fill='both', expand=True)
            ch['canvas'] = cv
            txt = ('Спайков: ' if is_ru else 'Spikes: ') + str(r['spikes'])
            txt += ' | ' + ('Частота: ' if is_ru else 'Firing rate: ')
            txt += '{} Hz'.format(r['firing_rate_hz'])
            if r['spikes'] == 0:
                txt += '  (' + ('ниже порога!' if is_ru else 'below threshold!') + ')'
            info_lbl.config(text=txt, fg='#1a7f37' if r['spikes'] > 0 else '#9a6700')

        tk.Button(dlg, text=('▶ Симулировать' if is_ru else '▶ Simulate'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run).pack(fill='x', padx=16, pady=(0, 8))
        run()

    def _open_msa(self):
        """Multiple sequence alignment viewer."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Множественное выравнивание' if is_ru else 'Multiple Alignment'))
        dlg.geometry('850x600')
        dlg.configure(bg=Theme.c('bg'))

        tk.Label(dlg, text=('Множественное выравнивание (MSA)'
                                       if is_ru else 'Multiple Sequence Alignment'),
                  bg=Theme.c('bg'), fg='#1a7f37',
                  font=('Arial', 15, 'bold')).pack(pady=10)
        tk.Label(dlg, text=('Вставьте последовательности (по одной на строку)'
                              if is_ru else 'Paste sequences (one per line)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack()

        inp = scrolledtext.ScrolledText(dlg, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                          font=('Courier', 10), height=6, relief='flat')
        inp.pack(fill='x', padx=16, pady=8)
        inp.insert('1.0', 'ATGGCAAGCTTAGCTAGCTAGCATCG\nATGGCAAGCTAAGCTAGCTAGCATCG\nATGGCTAGCTTAGCTAGGTAGCATCG\nATGGCAAGCTTAGCTACCTAGCATCG')

        out = scrolledtext.ScrolledText(dlg, bg=Theme.c('bg'), fg=Theme.c('text'),
                                         font=('Courier', 10), height=14, relief='flat',
                                         state='disabled')
        out.pack(fill='both', expand=True, padx=16, pady=8)

        def run():
            text = inp.get('1.0', 'end').strip()
            seqs = [s.strip() for s in text.split('\n') if s.strip()]
            if len(seqs) < 2:
                return
            r = msa_progressive(seqs)
            cons = msa_column_conservation(r['aligned'])
            out.config(state='normal')
            out.delete('1.0', 'end')
            for tag, color in [('match', '#1a7f37'), ('weak', '#9a6700'),
                                 ('low', '#57606a'), ('cons', '#0969da')]:
                out.tag_config(tag, foreground=color)
            for i, seq in enumerate(r['aligned']):
                out.insert('end', '{:8s} '.format(r['names'][i][:8]))
                for j, c in enumerate(seq):
                    score = cons[j] if j < len(cons) else 0
                    tag = 'match' if score >= 0.9 else 'weak' if score >= 0.6 else 'low'
                    out.insert('end', c, tag)
                out.insert('end', '\n')
            out.insert('end', '\n{:8s} '.format('CONS'))
            out.insert('end', r['consensus'], 'cons')
            out.insert('end', '\n\n')
            out.insert('end', ('Идентичность: ' if is_ru else 'Identity: ')
                        + '{}%\n'.format(r['identity']), 'cons')
            out.config(state='disabled')

        tk.Button(dlg, text=('▶ Выровнять' if is_ru else '▶ Align'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run).pack(fill='x', padx=16, pady=(0, 8))
        run()

    def _open_adme(self):
        """ADME / drug-likeness analyzer."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('ADME')
        dlg.geometry('600x600')
        dlg.configure(bg=Theme.c('bg'))

        tk.Label(dlg, text=('Анализ ADME / лекарственность'
                                       if is_ru else 'ADME / Drug-likeness'),
                  bg=Theme.c('bg'), fg='#8250df',
                  font=('Arial', 15, 'bold')).pack(pady=10)
        tk.Label(dlg, text=('Введите SMILES или параметры молекулы'
                              if is_ru else 'Enter SMILES or molecular parameters'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack()

        srow = tk.Frame(dlg, bg=Theme.c('bg')); srow.pack(fill='x', padx=16, pady=8)
        tk.Label(srow, text='SMILES:', bg=Theme.c('bg'), fg=Theme.c('text'),
                  font=('Arial', 9)).pack(side='left')
        smiles_var = tk.StringVar(value='CC(=O)Oc1ccccc1C(=O)O')
        tk.Entry(srow, textvariable=smiles_var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                  relief='flat').pack(side='left', fill='x', expand=True, padx=6)

        out = scrolledtext.ScrolledText(dlg, bg=Theme.c('bg'), fg=Theme.c('text'),
                                         font=('Arial', 10), height=20, relief='flat',
                                         state='disabled')
        out.pack(fill='both', expand=True, padx=16, pady=8)

        def run():
            smiles = smiles_var.get().strip()
            if not _ensure_rdkit():
                self._adme_show(out, None, 'RDKit не установлен' if is_ru else 'RDKit not installed')
                return
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    self._adme_show(out, None, 'Неверный SMILES' if is_ru else 'Invalid SMILES')
                    return
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                hd = Lipinski.NumHDonors(mol)
                ha = Lipinski.NumHAcceptors(mol)
                rb = Lipinski.NumRotatableBonds(mol)
                tpsa = Descriptors.TPSA(mol)
                r = adme_properties(mw, logp, hd, ha, rb, tpsa)
                r['_mw'] = mw; r['_logp'] = logp; r['_hd'] = hd
                r['_ha'] = ha; r['_rb'] = rb; r['_tpsa'] = tpsa
                self._adme_show(out, r, None)
            except Exception as e:
                biosim_log('warn', 'run: {}'.format(e))
                self._adme_show(out, None, str(e)[:80])

        tk.Button(dlg, text=('▶ Анализировать' if is_ru else '▶ Analyze'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run).pack(fill='x', padx=16, pady=(0, 8))
        run()

    def _adme_show(self, out, r, error):
        is_ru = (self.lang == 'ru')
        out.config(state='normal')
        out.delete('1.0', 'end')
        for tag, c in [('ok', '#1a7f37'), ('bad', '#cf222e'),
                        ('hdr', '#0969da'), ('dim', '#57606a')]:
            out.tag_config(tag, foreground=c, font=('Arial', 10, 'bold') if tag == 'hdr' else ('Arial', 10))
        if error:
            out.insert('end', '' + error + '\n', 'bad')
            out.config(state='disabled'); return
        out.insert('end', ('Свойства молекулы:\n' if is_ru else 'Molecular properties:\n'), 'hdr')
        out.insert('end', '  MW = {:.1f} g/mol\n'.format(r['_mw']))
        out.insert('end', '  LogP = {:.2f}\n'.format(r['_logp']))
        out.insert('end', '  H-доноры = {} | H-акцепторы = {}\n'.format(r['_hd'], r['_ha']))
        out.insert('end', '  Вращающиеся связи = {}\n'.format(r['_rb']))
        out.insert('end', '  TPSA = {:.1f} Å²\n\n'.format(r['_tpsa']))
        out.insert('end', ('Правило Липински:\n' if is_ru else 'Lipinski rule:\n'), 'hdr')
        for rule, ok in r['lipinski'].items():
            out.insert('end', '  {} {}\n'.format('' if ok else '', rule),
                        'ok' if ok else 'bad')
        out.insert('end', ('\nПравило Вебера:\n' if is_ru else '\nVeber rule:\n'), 'hdr')
        for rule, ok in r['veber'].items():
            out.insert('end', '  {} {}\n'.format('' if ok else '', rule),
                        'ok' if ok else 'bad')
        out.insert('end', '\n')
        out.insert('end', '  {} {}\n'.format(
            '' if r['druglike'] else '',
            ('Лекарственноподобна' if is_ru else 'Drug-like')),
            'ok' if r['druglike'] else 'bad')
        out.insert('end', '  {} {}\n'.format(
            '' if r['oral_bioavailable'] else '',
            ('Перорально биодоступна' if is_ru else 'Orally bioavailable')),
            'ok' if r['oral_bioavailable'] else 'bad')
        out.insert('end', '  {} {}\n'.format(
            '' if r['bbb_permeable'] else '',
            ('Проникает в мозг (ГЭБ)' if is_ru else 'BBB permeable')),
            'ok' if r['bbb_permeable'] else 'dim')
        out.config(state='disabled')

    def _open_epidemiology_pro(self):
        """Advanced epidemiology — interventions, age groups, cities."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Эпидемиология PRO' if is_ru else 'Epidemiology PRO'))
        dlg.geometry('950x700')
        dlg.configure(bg=Theme.c('bg'))

        tk.Label(dlg, text=('Продвинутая эпидемиология'
                                       if is_ru else 'Advanced Epidemiology'),
                  bg=Theme.c('bg'), fg='#cf222e',
                  font=('Arial', 15, 'bold')).pack(pady=8)

        nb = ttk.Notebook(dlg); nb.pack(fill='both', expand=True, padx=10, pady=6)

        # ─── TAB 1: Interventions (lockdown) ───
        tab1 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(tab1, text=('Локдаун' if is_ru else 'Lockdown'))
        tk.Label(tab1, text=('SEIR с вмешательствами (карантин, маски)'
                              if is_ru else 'SEIR with interventions'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        iv_ctrl = tk.Frame(tab1, bg=Theme.c('bg')); iv_ctrl.pack(fill='x', padx=12)
        iv_params = {}
        for i, (k, lbl, d) in enumerate([
                ('beta', 'β', 0.5), ('lockday', 'Локдаун с дня' if is_ru else 'Lockdown day', 30),
                ('lockend', 'до дня' if is_ru else 'until day', 70),
                ('strength', 'Сила % (сниж. β)' if is_ru else 'Strength %', 70)]):
            r = tk.Frame(iv_ctrl, bg=Theme.c('bg')); r.grid(row=i // 2, column=i % 2, padx=8, pady=3, sticky='w')
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9),
                      width=16, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=8, relief='flat').pack(side='left')
            iv_params[k] = v
        iv_fig = tk.Frame(tab1, bg=Theme.c('bg')); iv_fig.pack(fill='both', expand=True, padx=12, pady=6)
        iv_info = tk.Label(tab1, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        iv_info.pack()
        iv_ch = {'c': None}

        def run_iv():
            try:
                beta = float(iv_params['beta'].get())
                ld = float(iv_params['lockday'].get())
                le = float(iv_params['lockend'].get())
                strength = float(iv_params['strength'].get()) / 100.0
            except ValueError as _exc:
                biosim_log('warn', 'run_iv: {}'.format(_exc))
                return
            r = seir_with_interventions(10000, 10, beta, 0.2, 0.1, 150,
                interventions=[{'day': ld, 'end_day': le, 'beta_mult': 1 - strength}])
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if iv_ch['c']: iv_ch['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['I'], color='#cf222e', label='I', linewidth=1.5)
            ax.plot(r['times'], r['E'], color='#9a6700', label='E', alpha=0.7)
            ax.axvspan(ld, le, alpha=0.15, color='#0969da',
                        label=('Локдаун' if is_ru else 'Lockdown'))
            ax.set_xlabel(('Дни' if is_ru else 'Days'), color='#1c2128')
            ax.set_ylabel(('Заражённых' if is_ru else 'Infected'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=iv_fig); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); iv_ch['c'] = cv
            iv_info.config(text='R₀={} | {} {} {} {}'.format(
                r['R0'], ('пик' if is_ru else 'peak'), r['peak_infected'],
                ('на день' if is_ru else 'day'), r['peak_day']))
        tk.Button(tab1, text=('▶ Запустить' if is_ru else '▶ Run'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_iv).pack(fill='x', padx=12, pady=4)

        # ─── TAB 2: Cities (metapopulation) ───
        tab2 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(tab2, text=('Города' if is_ru else 'Cities'))
        tk.Label(tab2, text=('Распространение между городами'
                              if is_ru else 'Spread between cities'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        city_fig = tk.Frame(tab2, bg=Theme.c('bg')); city_fig.pack(fill='both', expand=True, padx=12, pady=6)
        city_info = tk.Label(tab2, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        city_info.pack()
        city_ch = {'c': None}

        def run_cities():
            cities = [{'name': 'A (старт)' if is_ru else 'A (start)', 'S': 10000, 'I': 20, 'R': 0},
                      {'name': 'B', 'S': 8000, 'I': 0, 'R': 0},
                      {'name': 'C', 'S': 5000, 'I': 0, 'R': 0},
                      {'name': 'D', 'S': 3000, 'I': 0, 'R': 0}]
            mig = [[0, 0.02, 0.01, 0.005], [0.02, 0, 0.015, 0.008],
                   [0.01, 0.015, 0, 0.02], [0.005, 0.008, 0.02, 0]]
            r = metapopulation_sir(cities, 0.4, 0.1, mig, 120)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if city_ch['c']: city_ch['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            colors = ['#cf222e', '#9a6700', '#1a7f37', '#0969da']
            for i, name in enumerate(r['names']):
                ax.plot(r['times'], r['cities'][name]['I'], color=colors[i % 4], label=name)
            ax.set_xlabel(('Дни' if is_ru else 'Days'), color='#1c2128')
            ax.set_ylabel(('Заражённых' if is_ru else 'Infected'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=city_fig); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); city_ch['c'] = cv
            peaks = ' | '.join('{}: {}д'.format(n, d) for n, d in r['peak_days'].items())
            city_info.config(text=('Пики: ' if is_ru else 'Peaks: ') + peaks)
        tk.Button(tab2, text=('▶ Запустить' if is_ru else '▶ Run'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_cities).pack(fill='x', padx=12, pady=4)

        # ─── TAB 3: Stochastic ───
        tab3 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(tab3, text=('Вспышки' if is_ru else 'Outbreaks'))
        tk.Label(tab3, text=('Стохастическая модель — вариабельность вспышек'
                              if is_ru else 'Stochastic model — outbreak variability'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        st_fig = tk.Frame(tab3, bg=Theme.c('bg')); st_fig.pack(fill='both', expand=True, padx=12, pady=6)
        st_info = tk.Label(tab3, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        st_info.pack()
        st_ch = {'c': None}

        def run_stoch():
            r = stochastic_sir(1000, 3, 0.25, 0.1, 80, n_runs=30)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if st_ch['c']: st_ch['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            for run in r['runs']:
                ax.plot(range(len(run)), run, color='#cf222e', alpha=0.25, linewidth=0.8)
            ax.set_xlabel(('Дни' if is_ru else 'Days'), color='#1c2128')
            ax.set_ylabel(('Заражённых' if is_ru else 'Infected'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=st_fig); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); st_ch['c'] = cv
            st_info.config(text='{} {} | {} {}'.format(
                ('Вероятность затухания:' if is_ru else 'Extinction prob:'),
                r['extinction_prob'],
                ('средний размер:' if is_ru else 'mean size:'), r['mean_final_size']))
        tk.Button(tab3, text=('▶ Запустить (30 прогонов)' if is_ru else '▶ Run (30 runs)'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_stoch).pack(fill='x', padx=12, pady=4)

        # ─── TAB 4: Age groups ───
        tab4 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(tab4, text=('Возрасты' if is_ru else 'Age groups'))
        tk.Label(tab4, text=('Эпидемия по возрастным группам (контактные матрицы)'
                              if is_ru else 'Epidemic by age groups (contact matrices)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        age_fig = tk.Frame(tab4, bg=Theme.c('bg')); age_fig.pack(fill='both', expand=True, padx=12, pady=6)
        age_info = tk.Label(tab4, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        age_info.pack()
        age_ch = {'c': None}

        def run_age():
            groups = [
                {'name': 'Дети' if is_ru else 'Children', 'S': 3000, 'I': 5, 'R': 0},
                {'name': 'Взрослые' if is_ru else 'Adults', 'S': 5000, 'I': 5, 'R': 0},
                {'name': 'Пожилые' if is_ru else 'Elderly', 'S': 2000, 'I': 0, 'R': 0},
            ]
            # Contact matrix: children contact each other most, elderly least
            cm = [[15, 5, 2], [5, 10, 3], [2, 3, 4]]
            r = age_structured_sir(groups, cm, 0.05, 0.1, 100)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if age_ch['c']: age_ch['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            colors = ['#1a7f37', '#0969da', '#cf222e']
            for i, name in enumerate(r['names']):
                ax.plot(r['times'], r['groups'][name]['I'], color=colors[i % 3], label=name)
            ax.set_xlabel(('Дни' if is_ru else 'Days'), color='#1c2128')
            ax.set_ylabel(('Заражённых' if is_ru else 'Infected'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=age_fig); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); age_ch['c'] = cv
            # Also show Rt example
            rt = effective_reproduction_number(4.0, 0.5)
            age_info.config(text='{} {}% | Rt({}) = {}'.format(
                ('Общий охват:' if is_ru else 'Attack rate:'), r['total_attack_rate'],
                ('50% восприимч.' if is_ru else '50% susceptible'), rt))
        tk.Button(tab4, text=('▶ Запустить' if is_ru else '▶ Run'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_age).pack(fill='x', padx=12, pady=4)

        run_iv()

    def _open_neuro_network(self):
        """Neural network simulator with spike raster."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Сеть нейронов' if is_ru else 'Neural Network'))
        dlg.geometry('900x680')
        dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Сеть нейронов (LIF) — растр спайков'
                                       if is_ru else 'Neural network (LIF) — spike raster'),
                  bg=Theme.c('bg'), fg='#8250df', font=('Arial', 15, 'bold')).pack(pady=8)

        ctrl = tk.Frame(dlg, bg=Theme.c('bg')); ctrl.pack(fill='x', padx=14, pady=6)
        params = {}
        for i, (k, lbl, d) in enumerate([
                ('n', 'Нейронов' if is_ru else 'Neurons', 30),
                ('conn', 'Связность %' if is_ru else 'Connectivity %', 20),
                ('cur', 'Ток I' if is_ru else 'Current I', 20)]):
            r = tk.Frame(ctrl, bg=Theme.c('bg')); r.pack(side='left', padx=10)
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=6, relief='flat').pack(side='left')
            params[k] = v
        fig_fr = tk.Frame(dlg, bg=Theme.c('bg')); fig_fr.pack(fill='both', expand=True, padx=14, pady=6)
        info = tk.Label(dlg, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info.pack()
        ch = {'c': None}

        def run():
            try:
                n = int(params['n'].get())
                conn = float(params['conn'].get()) / 100.0
                cur = float(params['cur'].get())
            except ValueError as _exc:
                biosim_log('warn', 'run: {}'.format(_exc))
                return
            r = lif_network(n, conn, I_ext=cur, t_max=200, seed=1)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch['c']: ch['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4.3), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            for ni, spikes in enumerate(r['spike_times']):
                ax.scatter(spikes, [ni] * len(spikes), s=4, color='#8250df')
            ax.set_xlabel(('Время (мс)' if is_ru else 'Time (ms)'), color='#1c2128')
            ax.set_ylabel(('Нейрон №' if is_ru else 'Neuron #'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.15)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=fig_fr); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch['c'] = cv
            info.config(text='{} {} | {} {} Hz'.format(
                ('Спайков:' if is_ru else 'Spikes:'), r['total_spikes'],
                ('средняя частота:' if is_ru else 'mean rate:'), r['mean_firing_rate_hz']))
        tk.Button(dlg, text=('▶ Симулировать' if is_ru else '▶ Simulate'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run).pack(fill='x', padx=14, pady=4)
        run()

    def _open_docking(self):
        """Molecular docking from SMILES."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Молекулярный докинг' if is_ru else 'Molecular Docking'))
        dlg.geometry('640x560')
        dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Молекулярный докинг'
                                       if is_ru else 'Molecular Docking'),
                  bg=Theme.c('bg'), fg='#0969da', font=('Arial', 15, 'bold')).pack(pady=8)
        tk.Label(dlg, text=('Связывание лиганда с белком-мишенью'
                              if is_ru else 'Ligand binding to target protein'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack()

        srow = tk.Frame(dlg, bg=Theme.c('bg')); srow.pack(fill='x', padx=16, pady=8)
        tk.Label(srow, text='SMILES:', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
        smiles_var = tk.StringVar(value='CC(=O)Oc1ccccc1C(=O)O')
        tk.Entry(srow, textvariable=smiles_var, bg=Theme.c('panel2'), fg=Theme.c('text'),
                  relief='flat').pack(side='left', fill='x', expand=True, padx=6)

        out = scrolledtext.ScrolledText(dlg, bg=Theme.c('bg'), fg=Theme.c('text'),
                                         font=('Courier', 10), height=18, relief='flat',
                                         state='disabled')
        out.pack(fill='both', expand=True, padx=16, pady=8)

        rot_var = tk.BooleanVar(value=True)
        tk.Checkbutton(dlg, text=('Перебирать ориентации лиганда (8 поворотов) — '
                                  'медленнее, но точнее' if is_ru else
                                  'Sample ligand orientations (8 rotations) — '
                                  'slower but more accurate'),
                       variable=rot_var, bg=Theme.c('bg'), fg=Theme.c('text'),
                       selectcolor=Theme.c('panel2'), font=('Arial', 9),
                       anchor='w').pack(fill='x', padx=16, pady=(0, 4))

        def run():
            out.config(state='normal'); out.delete('1.0', 'end')
            for tag, c in [('ok', '#1a7f37'), ('bad', '#cf222e'), ('hdr', '#0969da')]:
                out.tag_config(tag, foreground=c)
            if not _ensure_rdkit():
                out.insert('end', ('RDKit не установлен\n' if is_ru else 'RDKit not installed\n'), 'bad')
                out.config(state='disabled'); return
            out.insert('end', ('⏳ Генерация 3D и докинг...\n' if is_ru else '⏳ Generating 3D and docking...\n'))
            out.config(state='disabled'); out.update()
            r = dock_smiles(smiles_var.get().strip(),
                            sample_rotations=8 if rot_var.get() else 0)
            out.config(state='normal'); out.delete('1.0', 'end')
            for tag, c in [('ok', '#1a7f37'), ('bad', '#cf222e'), ('hdr', '#0969da')]:
                out.tag_config(tag, foreground=c)
            if 'error' in r:
                out.insert('end', '' + r['error'] + '\n', 'bad')
                out.config(state='disabled'); return
            out.insert('end', ('Результат докинга:\n' if is_ru else 'Docking result:\n'), 'hdr')
            out.insert('end', '  ' + ('Атомов лиганда: ' if is_ru else 'Ligand atoms: ') + str(r['ligand_atoms']) + '\n')
            out.insert('end', '  ' + ('Атомов рецептора: ' if is_ru else 'Receptor atoms: ') + str(r['receptor_atoms']) + '\n')
            out.insert('end', '  ' + ('Проверено поз: ' if is_ru else 'Poses tested: ') + str(r['poses_tested']) + '\n\n')
            score = r['best_score']
            out.insert('end', '  ' + ('Лучший score: ' if is_ru else 'Best score: ') + str(score) + '\n',
                        'ok' if score < 0 else 'bad')
            d = r.get('details', {})
            if d:
                out.insert('end', '  ' + ('Контактов: ' if is_ru else 'Contacts: ') + str(d.get('contacts', 0)) + '\n')
                out.insert('end', '  ' + ('Водородных связей: ' if is_ru else 'H-bonds: ') + str(d.get('h_bonds', 0)) + '\n')
                out.insert('end', '  ' + ('Качество: ' if is_ru else 'Quality: ') + str(d.get('binding_quality', '?')) + '\n', 'ok')
            if r.get('rotations_tried'):
                out.insert('end', '  ' + ('Ориентаций проверено: ' if is_ru
                                          else 'Orientations tried: ')
                           + str(r['rotations_tried'])
                           + ('  · лучшая при повороте {}°\n'.format(r.get('rotation_deg', 0))
                              if is_ru else
                              '  · best at {}°\n'.format(r.get('rotation_deg', 0))))
            kd = r.get('estimated_kd_nm', 0)
            if kd < 99999:
                out.insert('end', '\n  ' + ('Прибл. Kd: ' if is_ru else 'Est. Kd: ') + '{} nM\n'.format(kd))
            out.config(state='disabled')

        tk.Button(dlg, text=('▶ Докинг' if is_ru else '▶ Dock'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run).pack(fill='x', padx=16, pady=(0, 8))

    def _open_pk_advanced(self):
        """Advanced PK — 2-compartment and multiple dosing."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Фармакокинетика PRO' if is_ru else 'Advanced PK'))
        dlg.geometry('880x650'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Многокамерная фармакокинетика'
                                       if is_ru else 'Multi-compartment PK'),
                  bg=Theme.c('bg'), fg='#8250df', font=('Arial', 15, 'bold')).pack(pady=8)

        nb = ttk.Notebook(dlg); nb.pack(fill='both', expand=True, padx=10, pady=6)

        # Tab 1: 2-compartment
        t1 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t1, text=('2 камеры' if is_ru else '2-compartment'))
        p1 = {}
        c1 = tk.Frame(t1, bg=Theme.c('bg')); c1.pack(fill='x', padx=12, pady=6)
        for i, (k, lbl, d) in enumerate([('dose', 'Доза мг' if is_ru else 'Dose mg', 100),
                ('ka', 'ka (всасыв.)' if is_ru else 'ka', 1.0), ('k10', 'k10 (вывод)', 0.3),
                ('k12', 'k12', 0.2), ('k21', 'k21', 0.1), ('V1', 'V1 (л)', 10)]):
            r = tk.Frame(c1, bg=Theme.c('bg')); r.grid(row=i // 3, column=i % 3, padx=6, pady=3, sticky='w')
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9),
                      width=12, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=7, relief='flat').pack(side='left')
            p1[k] = v
        f1 = tk.Frame(t1, bg=Theme.c('bg')); f1.pack(fill='both', expand=True, padx=12, pady=6)
        ch1 = {'c': None}

        def run_2comp():
            try:
                dose = float(p1['dose'].get()); ka = float(p1['ka'].get())
                k10 = float(p1['k10'].get()); k12 = float(p1['k12'].get())
                k21 = float(p1['k21'].get()); V1 = float(p1['V1'].get())
            except ValueError as _exc:
                biosim_log('warn', 'run_2comp: {}'.format(_exc))
                return
            times = [i * 0.5 for i in range(49)]
            r = pk_two_compartment(dose, ka, k10, k12, k21, V1, times)
            if 'error' in r:
                return
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch1['c']: ch1['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 3.8), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['central'], color='#0969da', label=('Кровь' if is_ru else 'Central'))
            ax.plot(r['times'], r['peripheral'], color='#cf222e', label=('Ткани' if is_ru else 'Peripheral'))
            ax.set_xlabel(('Время (ч)' if is_ru else 'Time (h)'), color='#1c2128')
            ax.set_ylabel(('Концентрация' if is_ru else 'Concentration'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f1); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch1['c'] = cv
        tk.Button(t1, text=('▶ Рассчитать' if is_ru else '▶ Calculate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_2comp).pack(fill='x', padx=12, pady=4)

        # Tab 2: Multiple dosing
        t2 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t2, text=('Многократно' if is_ru else 'Multiple doses'))
        p2d = {}
        c2 = tk.Frame(t2, bg=Theme.c('bg')); c2.pack(fill='x', padx=12, pady=6)
        for i, (k, lbl, d) in enumerate([('dose', 'Доза мг' if is_ru else 'Dose', 100),
                ('ka', 'ka', 1.0), ('ke', 'ke', 0.15), ('V', 'V (л)', 10),
                ('tau', 'Интервал ч' if is_ru else 'Interval h', 8),
                ('n', 'Доз' if is_ru else 'Doses', 6)]):
            r = tk.Frame(c2, bg=Theme.c('bg')); r.grid(row=i // 3, column=i % 3, padx=6, pady=3, sticky='w')
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9),
                      width=12, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=7, relief='flat').pack(side='left')
            p2d[k] = v
        f2 = tk.Frame(t2, bg=Theme.c('bg')); f2.pack(fill='both', expand=True, padx=12, pady=6)
        info2 = tk.Label(t2, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info2.pack()
        ch2 = {'c': None}

        def run_multi():
            try:
                dose = float(p2d['dose'].get()); ka = float(p2d['ka'].get())
                ke = float(p2d['ke'].get()); V = float(p2d['V'].get())
                tau = float(p2d['tau'].get()); n = int(p2d['n'].get())
            except ValueError as _exc:
                biosim_log('warn', 'run_multi: {}'.format(_exc))
                return
            times = [i * 0.25 for i in range(int(n * tau / 0.25) + 20)]
            r = pk_multiple_dosing(dose, ka, ke, V, tau, n, times)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch2['c']: ch2['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 3.6), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['conc'], color='#8250df')
            ax.set_xlabel(('Время (ч)' if is_ru else 'Time (h)'), color='#1c2128')
            ax.set_ylabel(('Концентрация' if is_ru else 'Concentration'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f2); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch2['c'] = cv
            info2.config(text='Css_max={} | {} {}'.format(
                r['css_max'], ('накопление' if is_ru else 'accumulation'),
                r['accumulation_ratio']))
        tk.Button(t2, text=('▶ Рассчитать' if is_ru else '▶ Calculate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_multi).pack(fill='x', padx=12, pady=4)

        # ── 3 compartments ──────────────────────────────────────────
        # pk_three_compartment existed in the engine but nothing could reach
        # it, so drugs with deep-tissue distribution had no model here.
        t3 = tk.Frame(nb, bg=Theme.c('bg'))
        nb.add(t3, text=('3 камеры' if is_ru else '3-compartment'))
        tk.Label(t3, text=('Для препаратов с медленным депо — жир, кость, '
                           'мышца (многие анестетики)' if is_ru else
                           'For drugs with a slow deep compartment — fat, bone, '
                           'muscle (many anaesthetics)'),
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9),
                 wraplength=760, justify='left').pack(anchor='w', padx=12, pady=(6, 0))
        p3 = {}
        c3 = tk.Frame(t3, bg=Theme.c('bg')); c3.pack(fill='x', padx=12, pady=6)
        for i, (k, lbl, d) in enumerate([
                ('dose', 'Доза мг' if is_ru else 'Dose mg', 100),
                ('ka', 'ka (всасыв.)' if is_ru else 'ka', 1.0),
                ('k10', 'k10 (вывод)', 0.3), ('k12', 'k12', 0.2),
                ('k21', 'k21', 0.1), ('k13', 'k13 (глуб.)', 0.05),
                ('k31', 'k31', 0.01), ('V1', 'V1 (л)', 10)]):
            r = tk.Frame(c3, bg=Theme.c('bg'))
            r.grid(row=i // 4, column=i % 4, padx=6, pady=3, sticky='w')
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'),
                     font=('Arial', 9), width=12, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'),
                     width=7, relief='flat').pack(side='left')
            p3[k] = v
        f3 = tk.Frame(t3, bg=Theme.c('bg'))
        f3.pack(fill='both', expand=True, padx=12, pady=6)
        ch3 = {'c': None}
        info3 = tk.Label(t3, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                         font=('Arial', 10), wraplength=760, justify='left')
        info3.pack(fill='x', padx=12)

        def run_3comp():
            try:
                r = pk_three_compartment(
                    safe_float(p3['dose'].get(), 100), safe_float(p3['ka'].get(), 1.0),
                    safe_float(p3['k10'].get(), 0.3), safe_float(p3['k12'].get(), 0.2),
                    safe_float(p3['k21'].get(), 0.1), safe_float(p3['k13'].get(), 0.05),
                    safe_float(p3['k31'].get(), 0.01), safe_float(p3['V1'].get(), 10))
            except Exception as exc:
                biosim_log('error', 'pk_three_compartment: {}'.format(exc))
                info3.config(text=str(exc)); return
            try:
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                if ch3['c']: ch3['c'].get_tk_widget().destroy()
                fig = plt.Figure(figsize=(8, 3.6), facecolor='#f5f7fa')
                ax = fig.add_subplot(111, facecolor='#ffffff')
                ax.plot(r['times'], r['central'], color='#cf222e', lw=2,
                        label=('Кровь (центральная)' if is_ru else 'Central (blood)'))
                ax.plot(r['times'], r['peripheral1'], color='#0969da', lw=1.5, ls='--',
                        label=('Ткань 1 (быстрая)' if is_ru else 'Tissue 1 (fast)'))
                ax.plot(r['times'], r['peripheral2'], color='#1a7f37', lw=1.5, ls=':',
                        label=('Ткань 2 (глубокая)' if is_ru else 'Tissue 2 (deep)'))
                ax.set_xlabel(('Время (ч)' if is_ru else 'Time (h)'), color='#1c2128')
                ax.set_ylabel(('Концентрация' if is_ru else 'Concentration'),
                              color='#1c2128')
                ax.legend(fontsize=8); ax.grid(alpha=0.2)
                ax.tick_params(colors='#57606a')
                for s in ax.spines.values(): s.set_color('#d0d7de')
                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=f3); cv.draw()
                cv.get_tk_widget().pack(fill='both', expand=True); ch3['c'] = cv
            except Exception as exc:
                biosim_log('warn', 'run_3comp chart: {}'.format(exc))
            info3.config(text=(
                'Глубокая камера освобождается медленно — именно она задаёт '
                'длинный «хвост» после отмены препарата.' if is_ru else
                'The deep compartment empties slowly — it is what produces the '
                'long tail after dosing stops.'))

        tk.Button(t3, text=('▶ Рассчитать' if is_ru else '▶ Calculate'),
                  bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                  relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                  pady=4, command=run_3comp).pack(fill='x', padx=12, pady=4)

        # ── first-pass metabolism ───────────────────────────────────
        t4 = tk.Frame(nb, bg=Theme.c('bg'))
        nb.add(t4, text=('Первое прохождение' if is_ru else 'First pass'))
        tk.Label(t4, text=('Сколько препарата переживёт печень по дороге из '
                           'кишечника в кровь (модель хорошо перемешанной печени)'
                           if is_ru else
                           'How much drug survives the liver on the way from gut '
                           'to blood (well-stirred liver model)'),
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9),
                 wraplength=760, justify='left').pack(anchor='w', padx=12, pady=(6, 0))
        p4 = {}
        c4 = tk.Frame(t4, bg=Theme.c('bg')); c4.pack(fill='x', padx=12, pady=6)
        for i, (k, lbl, d) in enumerate([
                ('dose', 'Доза мг' if is_ru else 'Dose mg', 100),
                ('Q', 'Кровоток печени л/ч' if is_ru else 'Hepatic flow L/h', 90),
                ('CLint', 'Внутр. клиренс л/ч' if is_ru else 'Intrinsic CL L/h', 30)]):
            r = tk.Frame(c4, bg=Theme.c('bg'))
            r.grid(row=0, column=i, padx=8, pady=3, sticky='w')
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'),
                     font=('Arial', 9), width=20, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'),
                     width=8, relief='flat').pack(side='left')
            p4[k] = v
        f4 = tk.Frame(t4, bg=Theme.c('bg'))
        f4.pack(fill='both', expand=True, padx=12, pady=6)
        ch4 = {'c': None}
        info4 = tk.Label(t4, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                         font=('Arial', 11), wraplength=760, justify='left')
        info4.pack(fill='x', padx=12)

        def run_firstpass():
            dose = safe_float(p4['dose'].get(), 100)
            Q = safe_float(p4['Q'].get(), 90)
            try:
                r = first_pass_effect(dose, Q, safe_float(p4['CLint'].get(), 30))
            except Exception as exc:
                biosim_log('error', 'first_pass_effect: {}'.format(exc))
                info4.config(text=str(exc)); return
            e = r['extraction_ratio']
            verdict = (('Печень забирает {:.0f}% дозы. До крови доходит {} мг '
                        'из {:.0f}.').format(e * 100, r['systemic_dose'], dose)
                       if is_ru else
                       ('The liver removes {:.0f}% of the dose. {} mg of {:.0f} '
                        'reaches the blood.').format(e * 100, r['systemic_dose'], dose))
            if e > 0.7:
                verdict += (' Извлечение высокое — приём внутрь малоэффективен, '
                            'нужен другой путь введения.' if is_ru else
                            ' Extraction is high — the oral route is poor, another '
                            'route is needed.')
            cls = r.get('classification')
            if cls:
                verdict += ('  Извлечение: {}.'.format(cls) if is_ru
                            else '  Extraction: {}.'.format(cls))
            info4.config(text=verdict, fg='#cf222e' if e > 0.7 else '#1a7f37')
            # sweep CLint so the user sees where the cliff is
            try:
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                if ch4['c']: ch4['c'].get_tk_widget().destroy()
                xs = [Q * f for f in [0.02 * i for i in range(1, 101)]]
                ys = []
                for cl in xs:
                    ys.append(first_pass_effect(dose, Q, cl)['hepatic_availability'] * 100)
                fig = plt.Figure(figsize=(8, 3.2), facecolor='#f5f7fa')
                ax = fig.add_subplot(111, facecolor='#ffffff')
                ax.plot(xs, ys, color='#8250df', lw=2)
                ax.axvline(safe_float(p4['CLint'].get(), 30), color='#cf222e',
                           ls='--', lw=1)
                ax.set_xlabel(('Внутренний клиренс, л/ч' if is_ru
                               else 'Intrinsic clearance, L/h'), color='#1c2128')
                ax.set_ylabel(('Доходит до крови, %' if is_ru
                               else 'Reaches blood, %'), color='#1c2128')
                ax.grid(alpha=0.2); ax.tick_params(colors='#57606a')
                for s in ax.spines.values(): s.set_color('#d0d7de')
                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=f4); cv.draw()
                cv.get_tk_widget().pack(fill='both', expand=True); ch4['c'] = cv
            except Exception as exc:
                biosim_log('warn', 'run_firstpass chart: {}'.format(exc))

        tk.Button(t4, text=('▶ Рассчитать' if is_ru else '▶ Calculate'),
                  bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                  relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                  pady=4, command=run_firstpass).pack(fill='x', padx=12, pady=4)

        # ── bioequivalence ──────────────────────────────────────────
        t5 = tk.Frame(nb, bg=Theme.c('bg'))
        nb.add(t5, text=('Биоэквивалентность' if is_ru else 'Bioequivalence'))
        tk.Label(t5, text=('Сравнение дженерика с оригиналом по правилу FDA: '
                           'AUC и Cmax должны укладываться в 80–125%'
                           if is_ru else
                           'Generic versus reference by the FDA rule: AUC and '
                           'Cmax must fall within 80–125%'),
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9),
                 wraplength=760, justify='left').pack(anchor='w', padx=12, pady=(6, 0))
        p5 = {}
        c5 = tk.Frame(t5, bg=Theme.c('bg')); c5.pack(fill='x', padx=12, pady=6)
        for i, (k, lbl, d) in enumerate([
                ('dose', 'Доза мг' if is_ru else 'Dose mg', 100),
                ('ka_r', 'ka оригинала' if is_ru else 'ka reference', 1.0),
                ('ke_r', 'ke оригинала' if is_ru else 'ke reference', 0.25),
                ('ka_t', 'ka дженерика' if is_ru else 'ka test', 0.7),
                ('ke_t', 'ke дженерика' if is_ru else 'ke test', 0.25),
                ('V', 'V (л)', 20)]):
            r = tk.Frame(c5, bg=Theme.c('bg'))
            r.grid(row=i // 3, column=i % 3, padx=6, pady=3, sticky='w')
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'),
                     font=('Arial', 9), width=16, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'),
                     width=7, relief='flat').pack(side='left')
            p5[k] = v
        f5 = tk.Frame(t5, bg=Theme.c('bg'))
        f5.pack(fill='both', expand=True, padx=12, pady=6)
        ch5 = {'c': None}
        info5 = tk.Label(t5, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                         font=('Arial', 11), wraplength=760, justify='left')
        info5.pack(fill='x', padx=12)

        def run_be():
            dose = safe_float(p5['dose'].get(), 100)
            V = safe_float(p5['V'].get(), 20)
            times = [i * 0.25 for i in range(0, 97)]
            try:
                ref = pk_one_compartment(dose, safe_float(p5['ka_r'].get(), 1.0),
                                         safe_float(p5['ke_r'].get(), 0.25), V,
                                         times, route='oral')
                test = pk_one_compartment(dose, safe_float(p5['ka_t'].get(), 0.7),
                                          safe_float(p5['ke_t'].get(), 0.25), V,
                                          times, route='oral')
                ref = list(ref); test = list(test)
                r = bioequivalence(times, test, times, ref)
            except Exception as exc:
                biosim_log('error', 'bioequivalence: {}'.format(exc))
                info5.config(text=str(exc)); return
            ok = r.get('bioequivalent', False)
            info5.config(
                text=(('ЭКВИВАЛЕНТНЫ' if ok else 'НЕ эквивалентны') if is_ru
                      else ('BIOEQUIVALENT' if ok else 'NOT bioequivalent'))
                + ' · AUC {:.0f}% · Cmax {:.0f}%'.format(
                    r.get('auc_ratio', 0), r.get('cmax_ratio', 0)),
                fg='#1a7f37' if ok else '#cf222e')
            try:
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                if ch5['c']: ch5['c'].get_tk_widget().destroy()
                fig = plt.Figure(figsize=(8, 3.2), facecolor='#f5f7fa')
                ax = fig.add_subplot(111, facecolor='#ffffff')
                ax.plot(times, ref, color='#0969da', lw=2,
                        label=('Оригинал' if is_ru else 'Reference'))
                ax.plot(times, test, color='#bc4c00', lw=2, ls='--',
                        label=('Дженерик' if is_ru else 'Test'))
                ax.set_xlabel(('Время (ч)' if is_ru else 'Time (h)'), color='#1c2128')
                ax.set_ylabel(('Концентрация' if is_ru else 'Concentration'),
                              color='#1c2128')
                ax.legend(fontsize=8); ax.grid(alpha=0.2)
                ax.tick_params(colors='#57606a')
                for s in ax.spines.values(): s.set_color('#d0d7de')
                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=f5); cv.draw()
                cv.get_tk_widget().pack(fill='both', expand=True); ch5['c'] = cv
            except Exception as exc:
                biosim_log('warn', 'run_be chart: {}'.format(exc))

        tk.Button(t5, text=('▶ Сравнить' if is_ru else '▶ Compare'),
                  bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                  relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                  pady=4, command=run_be).pack(fill='x', padx=12, pady=4)

        # ── genetic dose-response ───────────────────────────────────
        t6 = tk.Frame(nb, bg=Theme.c('bg'))
        nb.add(t6, text=('Доза-ответ (промотор)' if is_ru
                         else 'Dose-response (promoter)'))
        tk.Label(t6, text=('Сколько белка даст индуцируемый промотор при разной '
                           'концентрации индуктора — например IPTG. Коэффициент '
                           'Хилла задаёт, насколько резким будет переключение.'
                           if is_ru else
                           'How much protein an inducible promoter yields across '
                           'inducer concentrations — IPTG, for instance. The Hill '
                           'coefficient sets how sharp the switch is.'),
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9),
                 wraplength=760, justify='left').pack(anchor='w', padx=12, pady=(6, 0))
        p6 = {}
        c6 = tk.Frame(t6, bg=Theme.c('bg')); c6.pack(fill='x', padx=12, pady=6)
        for i, (k, lbl, d) in enumerate([
                ('vmax', 'Максимум', 100), ('km', 'Km (EC50)', 10),
                ('n', 'Коэф. Хилла' if is_ru else 'Hill coeff.', 2),
                ('leaky', 'Течь (базовый)' if is_ru else 'Leaky basal', 5)]):
            r = tk.Frame(c6, bg=Theme.c('bg'))
            r.grid(row=0, column=i, padx=8, pady=3, sticky='w')
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'),
                     font=('Arial', 9), width=14, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'),
                     width=7, relief='flat').pack(side='left')
            p6[k] = v
        f6 = tk.Frame(t6, bg=Theme.c('bg'))
        f6.pack(fill='both', expand=True, padx=12, pady=6)
        ch6 = {'c': None}
        info6 = tk.Label(t6, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                         font=('Arial', 10), wraplength=760, justify='left')
        info6.pack(fill='x', padx=12)

        def run_dr():
            km = safe_float(p6['km'].get(), 10)
            rng = [km * (10 ** (-2 + 0.05 * i)) for i in range(81)]
            try:
                r = dose_response_genetic(rng, safe_float(p6['vmax'].get(), 100),
                                          km, safe_float(p6['n'].get(), 2),
                                          safe_float(p6['leaky'].get(), 5))
            except Exception as exc:
                biosim_log('error', 'dose_response_genetic: {}'.format(exc))
                info6.config(text=str(exc)); return
            outs = r.get('output', [])
            try:
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                if ch6['c']: ch6['c'].get_tk_widget().destroy()
                fig = plt.Figure(figsize=(8, 3.4), facecolor='#f5f7fa')
                ax = fig.add_subplot(111, facecolor='#ffffff')
                ax.semilogx(rng, outs, color='#1a7f37', lw=2)
                ax.axvline(km, color='#8250df', ls='--', lw=1)
                ax.text(km, max(outs) * 0.1, ' EC50', color='#8250df', fontsize=8)
                ax.axhline(safe_float(p6['leaky'].get(), 5), color='#cf222e',
                           ls=':', lw=1)
                ax.set_xlabel(('Индуктор (усл. ед., лог. шкала)' if is_ru
                               else 'Inducer (a.u., log scale)'), color='#1c2128')
                ax.set_ylabel(('Выход белка' if is_ru else 'Protein output'),
                              color='#1c2128')
                ax.grid(alpha=0.2, which='both'); ax.tick_params(colors='#57606a')
                for s in ax.spines.values(): s.set_color('#d0d7de')
                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=f6); cv.draw()
                cv.get_tk_widget().pack(fill='both', expand=True); ch6['c'] = cv
            except Exception as exc:
                biosim_log('warn', 'run_dr chart: {}'.format(exc))
            n_val = safe_float(p6['n'].get(), 2)
            sharp = ('переключение резкое, почти как тумблер' if n_val >= 3 else
                     'переключение плавное — промежуточные уровни достижимы')
            fold = r.get('fold_change')
            info6.config(text=('EC50 = {} · размах включения ×{} · {}'.format(
                r.get('ec50', km), fold if fold else '—', sharp) if is_ru else
                'EC50 = {} · fold change x{}'.format(r.get('ec50', km),
                                                     fold if fold else '-')))

        tk.Button(t6, text=('▶ Построить' if is_ru else '▶ Plot'),
                  bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                  relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                  pady=4, command=run_dr).pack(fill='x', padx=12, pady=4)

        run_2comp()

    def _open_dose_calculator(self):
        """Dose calculator — by weight and renal function."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Калькулятор доз' if is_ru else 'Dose Calculator'))
        dlg.geometry('560x520'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Калькулятор дозировки' if is_ru else 'Dose Calculator'),
                  bg=Theme.c('bg'), fg='#8250df', font=('Arial', 15, 'bold')).pack(pady=8)

        form = tk.Frame(dlg, bg=Theme.c('bg')); form.pack(fill='x', padx=20, pady=10)
        vars = {}
        for k, lbl, d in [('dpk', 'Доза на кг (мг/кг)' if is_ru else 'Dose per kg', 5),
                          ('weight', 'Вес (кг)' if is_ru else 'Weight (kg)', 70),
                          ('crcl', 'Клиренс креатинина (мл/мин)' if is_ru else 'CrCl (mL/min)', 90)]:
            r = tk.Frame(form, bg=Theme.c('bg')); r.pack(fill='x', pady=5)
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 10),
                      width=28, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=10, relief='flat').pack(side='left')
            vars[k] = v

        out = scrolledtext.ScrolledText(dlg, bg=Theme.c('bg'), fg=Theme.c('text'),
                                         font=('Arial', 11), height=10, relief='flat', state='disabled')
        out.pack(fill='both', expand=True, padx=20, pady=10)

        def calc():
            try:
                dpk = float(vars['dpk'].get()); weight = float(vars['weight'].get())
                crcl = float(vars['crcl'].get())
            except ValueError as _exc:
                biosim_log('warn', 'calc: {}'.format(_exc))
                return
            base = dose_by_weight(dpk, weight)
            renal = renal_dose_adjustment(base, crcl)
            out.config(state='normal'); out.delete('1.0', 'end')
            out.tag_config('hdr', foreground='#0969da', font=('Arial', 11, 'bold'))
            out.tag_config('val', foreground='#1a7f37', font=('Arial', 12, 'bold'))
            out.insert('end', ('Доза по весу:\n' if is_ru else 'Weight-based dose:\n'), 'hdr')
            out.insert('end', '  {} мг ({} мг/кг × {} кг)\n\n'.format(base, dpk, weight), 'val')
            out.insert('end', ('С учётом почек:\n' if is_ru else 'Renal-adjusted:\n'), 'hdr')
            out.insert('end', '  {} мг\n'.format(renal['adjusted_dose']), 'val')
            out.insert('end', '  ' + ('функция почек: ' if is_ru else 'renal function: ')
                        + renal['renal_function'] + '\n')
            out.insert('end', '  ' + ('коэффициент: ' if is_ru else 'factor: ')
                        + str(renal['factor']) + '\n')
            out.config(state='disabled')
        tk.Button(dlg, text=('▶ Рассчитать' if is_ru else '▶ Calculate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=5, command=calc).pack(fill='x', padx=20, pady=(0, 10))
        calc()


    def _open_synbio_tools(self):
        """SynBio advanced tools — Gillespie, logic gates, codon optimize."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('SynBio инструменты' if is_ru else 'SynBio Tools'))
        dlg.geometry('820x650'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Инструменты синтетической биологии'
                                       if is_ru else 'Synthetic Biology Tools'),
                  bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 15, 'bold')).pack(pady=8)

        nb = ttk.Notebook(dlg); nb.pack(fill='both', expand=True, padx=10, pady=6)

        # Tab 1: Logic gates
        t1 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t1, text=('Логические гейты' if is_ru else 'Logic gates'))
        tk.Label(t1, text=('Генетические логические элементы (ДНК-компьютер)'
                              if is_ru else 'Genetic logic gates (DNA computing)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        gate_var = tk.StringVar(value='AND')
        grow = tk.Frame(t1, bg=Theme.c('bg')); grow.pack(pady=6)
        tk.Label(grow, text=('Гейт:' if is_ru else 'Gate:'), bg=Theme.c('bg'),
                  fg=Theme.c('text'), font=('Arial', 10)).pack(side='left')
        for g in ['AND', 'OR', 'NOT', 'NAND', 'NOR', 'XOR']:
            tk.Radiobutton(grow, text=g, variable=gate_var, value=g, bg=Theme.c('bg'),
                            fg='#0969da', selectcolor='#eef1f5', font=('Arial', 9)).pack(side='left', padx=3)
        gout = scrolledtext.ScrolledText(t1, bg=Theme.c('bg'), fg=Theme.c('text'),
                                          font=('Courier', 11), height=12, relief='flat', state='disabled')
        gout.pack(fill='both', expand=True, padx=12, pady=8)

        def run_gate():
            gt = gate_var.get()
            gout.config(state='normal'); gout.delete('1.0', 'end')
            gout.tag_config('hdr', foreground='#0969da', font=('Courier', 11, 'bold'))
            gout.tag_config('on', foreground='#1a7f37'); gout.tag_config('off', foreground=Theme.c('text_dim'))
            gout.insert('end', ('Таблица истинности ' if is_ru else 'Truth table ') + gt + ':\n\n', 'hdr')
            gout.insert('end', '  A  B  | OUT  ' + ('экспрессия' if is_ru else 'expression') + '\n')
            gout.insert('end', '  ─────────────────────\n')
            inputs = [(0, 0), (0, 1), (1, 0), (1, 1)] if gt not in ('NOT',) else [(0, None), (1, None)]
            for a, b in inputs:
                r = logic_gate_dna(gt, a, b)
                b_str = str(b) if b is not None else '-'
                tag = 'on' if r['output'] else 'off'
                gout.insert('end', '  {}  {}  |  {}   {:.2f}\n'.format(a, b_str, r['output'], r['expression']), tag)
            gout.config(state='disabled')
        tk.Button(t1, text=('▶ Показать' if is_ru else '▶ Show'), bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run_gate).pack(fill='x', padx=12, pady=4)

        # Tab 2: Codon optimization
        t2 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t2, text=('Оптимизация кодонов' if is_ru else 'Codon optimize'))
        tk.Label(t2, text=('Оптимизация ДНК для экспрессии в E. coli'
                              if is_ru else 'Optimize DNA for E. coli expression'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        tk.Label(t2, text=('Белок (одно-буквенный код):' if is_ru else 'Protein (1-letter):'),
                  bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(anchor='w', padx=12)
        codon_in = tk.Entry(t2, bg=Theme.c('panel2'), fg=Theme.c('text'), relief='flat', font=('Courier', 11))
        codon_in.pack(fill='x', padx=12, pady=4); codon_in.insert(0, 'MAKVLILACG')
        codon_out = scrolledtext.ScrolledText(t2, bg=Theme.c('bg'), fg=Theme.c('text'),
                                               font=('Courier', 10), height=10, relief='flat', state='disabled')
        codon_out.pack(fill='both', expand=True, padx=12, pady=8)

        def run_codon():
            protein = codon_in.get().strip().upper()
            if not protein:
                return
            r = codon_optimize_detailed(protein, is_protein=True)
            codon_out.config(state='normal'); codon_out.delete('1.0', 'end')
            codon_out.tag_config('hdr', foreground='#0969da', font=('Courier', 10, 'bold'))
            codon_out.tag_config('dna', foreground='#1a7f37')
            codon_out.insert('end', ('Оптимизированная ДНК:\n' if is_ru else 'Optimized DNA:\n'), 'hdr')
            codon_out.insert('end', r['optimized_dna'] + '\n\n', 'dna')
            codon_out.insert('end', ('Длина: ' if is_ru else 'Length: ') + str(r['length']) + ' bp\n')
            codon_out.insert('end', ('Белок: ' if is_ru else 'Protein: ') + r['protein'] + '\n')
            codon_out.config(state='disabled')
        tk.Button(t2, text=('▶ Оптимизировать' if is_ru else '▶ Optimize'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_codon).pack(fill='x', padx=12, pady=4)

        # Tab 3: Gillespie stochastic
        t3 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t3, text=('Стохастика (Gillespie)' if is_ru else 'Gillespie SSA'))
        tk.Label(t3, text=('Точная стохастическая симуляция реакций\n'
                              'Пример: A → B → C'
                              if is_ru else 'Exact stochastic simulation\nExample: A → B → C'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        f3 = tk.Frame(t3, bg=Theme.c('bg')); f3.pack(fill='both', expand=True, padx=12, pady=6)
        ch3 = {'c': None}

        def run_gillespie():
            reactions = [
                {'rate': 0.1, 'reactants': {'A': 1}, 'products': {'B': 1}},
                {'rate': 0.05, 'reactants': {'B': 1}, 'products': {'C': 1}},
            ]
            initial = {'A': 100, 'B': 0, 'C': 0}
            r = gillespie_ssa(reactions, initial, t_max=50, seed=1)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch3['c']: ch3['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            for sp, color in [('A', '#0969da'), ('B', '#9a6700'), ('C', '#1a7f37')]:
                if sp in r['species']:
                    ax.plot(r['times'], r['species'][sp], color=color, label=sp, drawstyle='steps-post')
            ax.set_xlabel(('Время' if is_ru else 'Time'), color='#1c2128')
            ax.set_ylabel(('Молекул' if is_ru else 'Molecules'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f3); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch3['c'] = cv
        tk.Button(t3, text=('▶ Симулировать' if is_ru else '▶ Simulate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_gillespie).pack(fill='x', padx=12, pady=4)

        run_gate()


    def _open_neuro_plasticity(self):
        """Synaptic plasticity (STDP) + FitzHugh-Nagumo neuron."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Пластичность и динамика' if is_ru else 'Plasticity & Dynamics'))
        dlg.geometry('860x650'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Синаптическая пластичность и нейродинамика'
                                       if is_ru else 'Synaptic plasticity & neurodynamics'),
                  bg=Theme.c('bg'), fg='#8250df', font=('Arial', 15, 'bold')).pack(pady=8)

        nb = ttk.Notebook(dlg); nb.pack(fill='both', expand=True, padx=10, pady=6)

        # Tab 1: STDP
        t1 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t1, text='STDP')
        tk.Label(t1, text=('Пластичность по времени спайков (обучение синапса)'
                              if is_ru else 'Spike-timing-dependent plasticity (synapse learning)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        f1 = tk.Frame(t1, bg=Theme.c('bg')); f1.pack(fill='both', expand=True, padx=12, pady=6)
        info1 = tk.Label(t1, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info1.pack()
        ch1 = {'c': None}

        def run_stdp():
            # Pre fires before post (potentiation)
            pre = [10, 30, 50, 70]
            post = [12, 32, 52, 72]
            r = synapse_stdp(pre, post)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch1['c']: ch1['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 3.8), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(range(len(r['w_history'])), r['w_history'], color='#8250df', marker='o', markersize=3)
            ax.set_xlabel(('Событие спайка' if is_ru else 'Spike event'), color='#1c2128')
            ax.set_ylabel(('Вес синапса' if is_ru else 'Synapse weight'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f1); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch1['c'] = cv
            info1.config(text='{} {} ({})'.format(
                ('Финальный вес:' if is_ru else 'Final weight:'),
                r['final_weight'],
                ('усилен' if r['direction'] == 'strengthened' else 'ослаблен') if is_ru else r['direction']))
        tk.Button(t1, text=('▶ Симулировать (pre→post)' if is_ru else '▶ Simulate (pre→post)'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_stdp).pack(fill='x', padx=12, pady=4)

        # Tab 2: FitzHugh-Nagumo
        t2 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t2, text='FitzHugh-Nagumo')
        tk.Label(t2, text=('Упрощённая модель нейрона (быстрее Ходжкина-Хаксли)'
                              if is_ru else 'Simplified neuron model (faster than HH)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        crow = tk.Frame(t2, bg=Theme.c('bg')); crow.pack(fill='x', padx=12, pady=4)
        tk.Label(crow, text=('Ток I:' if is_ru else 'Current I:'), bg=Theme.c('bg'),
                  fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
        fn_cur = tk.DoubleVar(value=0.5)
        tk.Scale(crow, from_=0, to=2, resolution=0.1, orient='horizontal', variable=fn_cur,
                  bg=Theme.c('bg'), fg=Theme.c('text'), highlightthickness=0, troughcolor='#eef1f5',
                  length=300).pack(side='left', padx=8)
        f2 = tk.Frame(t2, bg=Theme.c('bg')); f2.pack(fill='both', expand=True, padx=12, pady=6)
        ch2 = {'c': None}

        def run_fn():
            r = fitzhugh_nagumo(fn_cur.get(), t_max=100)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch2['c']: ch2['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 3.8), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['v'], color='#8250df', label='v (' + ('мембрана' if is_ru else 'membrane') + ')')
            ax.plot(r['times'], r['w'], color='#0969da', label='w (' + ('восстан.' if is_ru else 'recovery') + ')', alpha=0.7)
            ax.set_xlabel(('Время' if is_ru else 'Time'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f2); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch2['c'] = cv
        tk.Button(t2, text=('▶ Симулировать' if is_ru else '▶ Simulate'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_fn).pack(fill='x', padx=12, pady=4)

        run_stdp()

    def _open_phylo_builder(self):
        """Build phylogenetic tree from sequences (NJ + JC distance)."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Дерево NJ' if is_ru else 'NJ Tree'))
        dlg.geometry('760x620'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Филогенетическое дерево (Neighbor-Joining)'
                                       if is_ru else 'Phylogenetic tree (Neighbor-Joining)'),
                  bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 14, 'bold')).pack(pady=8)
        tk.Label(dlg, text=('Вставьте выровненные последовательности (по одной на строку)'
                              if is_ru else 'Paste aligned sequences (one per line)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack()

        inp = scrolledtext.ScrolledText(dlg, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                         font=('Courier', 10), height=6, relief='flat')
        inp.pack(fill='x', padx=16, pady=8)
        inp.insert('1.0', 'ATGGCAAGCTTAGCTAGCTAGC\nATGGCAAGCTAAGCTAGCTAGC\n'
                            'ATGGCTAGCTTAGCTAGGTAGC\nATGGCAAGCTTAGCTACCTAGC')

        out = scrolledtext.ScrolledText(dlg, bg=Theme.c('bg'), fg=Theme.c('text'),
                                         font=('Courier', 10), height=14, relief='flat', state='disabled')
        out.pack(fill='both', expand=True, padx=16, pady=8)

        def build():
            text = inp.get('1.0', 'end').strip()
            seqs = [s.strip() for s in text.split('\n') if s.strip()]
            if len(seqs) < 3:
                return
            names = ['Seq{}'.format(i + 1) for i in range(len(seqs))]
            D = build_distance_matrix(seqs)
            tree = neighbor_joining_tree(D, names)
            out.config(state='normal'); out.delete('1.0', 'end')
            out.tag_config('hdr', foreground='#0969da', font=('Courier', 10, 'bold'))
            out.tag_config('tree', foreground='#1a7f37')
            out.insert('end', ('Матрица расстояний (Jukes-Cantor):\n' if is_ru else 'Distance matrix (JC):\n'), 'hdr')
            out.insert('end', '       ' + '  '.join('{:>6}'.format(n) for n in names) + '\n')
            for i, name in enumerate(names):
                out.insert('end', '{:>6} '.format(name))
                for j in range(len(names)):
                    out.insert('end', '{:>7.3f} '.format(D[i][j]))
                out.insert('end', '\n')
            out.insert('end', '\n' + ('Дерево (Newick):\n' if is_ru else 'Tree (Newick):\n'), 'hdr')
            out.insert('end', tree['newick'] + '\n\n', 'tree')
            out.insert('end', ('Шаги объединения:\n' if is_ru else 'Merge steps:\n'), 'hdr')
            for a, b, da, db in tree['steps']:
                out.insert('end', '  {} + {} (ветви {:.3f}, {:.3f})\n'.format(
                    a[:15], b[:15], da, db))
            out.config(state='disabled')

        tk.Button(dlg, text=('▶ Построить дерево' if is_ru else '▶ Build tree'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=build).pack(fill='x', padx=16, pady=(0, 8))
        build()


    def _open_molbio_pro(self):
        """Professional MolBio tools — BLAST, plasmid annotation, RNA, primers."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Молекулярная биология PRO' if is_ru else 'MolBio PRO'))
        dlg.geometry('960x720'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Профессиональная молекулярная биология'
                                       if is_ru else 'Professional Molecular Biology'),
                  bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 15, 'bold')).pack(pady=8)

        nb = ttk.Notebook(dlg); nb.pack(fill='both', expand=True, padx=10, pady=6)

        # ═══ TAB 1: BLAST ═══
        t1 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t1, text='BLAST')
        tk.Label(t1, text=('Локальное выравнивание с e-value и bit-score'
                              if is_ru else 'Local alignment with e-value and bit-score'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        tk.Label(t1, text=('Запрос (query):' if is_ru else 'Query:'),
                  bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(anchor='w', padx=12)
        blast_q = tk.Entry(t1, bg=Theme.c('panel2'), fg=Theme.c('text'), relief='flat', font=('Courier', 10))
        blast_q.pack(fill='x', padx=12, pady=2); blast_q.insert(0, 'ATGGCAAGCTTAGCTAGC')
        tk.Label(t1, text=('База (по одной на строку):' if is_ru else 'Database (one per line):'),
                  bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(anchor='w', padx=12)
        blast_db = scrolledtext.ScrolledText(t1, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                              font=('Courier', 9), height=4, relief='flat')
        blast_db.pack(fill='x', padx=12, pady=2)
        blast_db.insert('1.0', 'ATGGCAAGCTTAGCTAGCTAGCATCG\nTTTGGGCCCAAATTTGGG\nATGGCAAGCTTACCCGGGTAG')
        blast_out = scrolledtext.ScrolledText(t1, bg=Theme.c('bg'), fg=Theme.c('text'),
                                              font=('Courier', 9), height=12, relief='flat', state='disabled')
        blast_out.pack(fill='both', expand=True, padx=12, pady=6)

        def run_blast():
            query = blast_q.get().strip()
            db = [s.strip() for s in blast_db.get('1.0', 'end').split('\n') if s.strip()]
            if not query or not db:
                return
            hits = blast_search(query, db)
            o = blast_out; o.config(state='normal'); o.delete('1.0', 'end')
            o.tag_config('hdr', foreground='#0969da', font=('Courier', 9, 'bold'))
            o.tag_config('match', foreground='#1a7f37')
            if not hits:
                o.insert('end', ('Совпадений не найдено\n' if is_ru else 'No hits found\n'))
            for h in hits:
                o.insert('end', '> {} | '.format(h['name']), 'hdr')
                o.insert('end', 'score={} bits={} E={:.2e} id={}%\n'.format(
                    h['score'], h['bit_score'], h['e_value'], h['identity']))
                o.insert('end', 'Q: {}\n'.format(h['q_aln'][:60]))
                o.insert('end', '   {}\n'.format(h['midline'][:60]), 'match')
                o.insert('end', 'S: {}\n\n'.format(h['s_aln'][:60]))
            o.config(state='disabled')
        tk.Button(t1, text=('▶ BLAST', None)[0] if False else '▶ BLAST', bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_blast).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 2: Plasmid Annotation ═══
        t2 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t2, text=('Аннотация' if is_ru else 'Annotation'))
        tk.Label(t2, text=('Автоматическая аннотация плазмиды (гены, промоторы, ori)'
                              if is_ru else 'Auto-annotate plasmid (genes, promoters, ori)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        plasmid_in = scrolledtext.ScrolledText(t2, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                                font=('Courier', 9), height=5, relief='flat')
        plasmid_in.pack(fill='x', padx=12, pady=4)
        plasmid_in.insert('1.0', 'TAATACGACTCACTATAGGGAGGAGGATG' + 'GCAGCT' * 60 + 'TAA')
        plasmid_out = scrolledtext.ScrolledText(t2, bg=Theme.c('bg'), fg=Theme.c('text'),
                                                 font=('Courier', 9), height=13, relief='flat', state='disabled')
        plasmid_out.pack(fill='both', expand=True, padx=12, pady=6)

        def run_annotate():
            seq = plasmid_in.get('1.0', 'end').strip()
            if len(seq) < 30:
                return
            ann = annotate_plasmid(seq, min_orf_len=60)
            o = plasmid_out; o.config(state='normal'); o.delete('1.0', 'end')
            o.tag_config('hdr', foreground='#0969da', font=('Courier', 9, 'bold'))
            colors = {'promoter': '#9a6700', 'ORF': '#1a7f37', 'RBS': '#0550ae',
                      'ori': '#cf222e', 'resistance': '#8250df'}
            for t, c in colors.items():
                o.tag_config(t, foreground=c)
            o.insert('end', 'Length: {} bp | GC: {}% | {} features\n\n'.format(
                ann['length'], ann['gc_content'], ann['feature_count']), 'hdr')
            for f in ann['features']:
                o.insert('end', '{:>5}-{:<5} '.format(f['start'], f['end']))
                o.insert('end', '[{}] '.format(f['strand']))
                o.insert('end', '{}: {}\n'.format(f['type'], f['name']),
                          f['type'] if f['type'] in colors else 'hdr')
            o.config(state='disabled')
        tk.Button(t2, text=('▶ Аннотировать' if is_ru else '▶ Annotate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_annotate).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 3: RNA Structure ═══
        t3 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t3, text=('РНК структура' if is_ru else 'RNA fold'))
        tk.Label(t3, text=('Вторичная структура РНК (алгоритм Нуссинова)'
                              if is_ru else 'RNA secondary structure (Nussinov)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        rna_in = tk.Entry(t3, bg=Theme.c('panel2'), fg=Theme.c('text'), relief='flat', font=('Courier', 11))
        rna_in.pack(fill='x', padx=12, pady=4)
        rna_in.insert(0, 'GGGAAACUCCGUUUCCC')
        rna_out = scrolledtext.ScrolledText(t3, bg=Theme.c('bg'), fg=Theme.c('text'),
                                            font=('Courier', 12), height=12, relief='flat', state='disabled')
        rna_out.pack(fill='both', expand=True, padx=12, pady=6)

        def run_rna():
            seq = rna_in.get().strip()
            if len(seq) < 4:
                return
            r = rna_secondary_structure(seq)
            o = rna_out; o.config(state='normal'); o.delete('1.0', 'end')
            o.tag_config('hdr', foreground='#0969da', font=('Courier', 11, 'bold'))
            o.tag_config('seq', foreground='#1a7f37', font=('Courier', 13))
            o.tag_config('struct', foreground='#9a6700', font=('Courier', 13))
            o.insert('end', ('Последовательность:\n' if is_ru else 'Sequence:\n'), 'hdr')
            o.insert('end', r['sequence'] + '\n', 'seq')
            o.insert('end', r['structure'] + '\n\n', 'struct')
            o.insert('end', ('Пар оснований: ' if is_ru else 'Base pairs: ')
                      + str(r['num_pairs']) + '\n')
            o.insert('end', ('Оценка MFE: ' if is_ru else 'MFE estimate: ')
                      + str(r['mfe_estimate']) + ' kcal/mol\n')
            o.config(state='disabled')
        tk.Button(t3, text=('▶ Свернуть' if is_ru else '▶ Fold'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_rna).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 4: Pro Primers ═══
        t4 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t4, text=('Праймеры PRO' if is_ru else 'Pro Primers'))
        tk.Label(t4, text=('Дизайн праймеров с точным Tm (nearest-neighbor) + QC'
                              if is_ru else 'Primer design with accurate NN Tm + QC'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        primer_in = scrolledtext.ScrolledText(t4, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                               font=('Courier', 9), height=4, relief='flat')
        primer_in.pack(fill='x', padx=12, pady=4)
        primer_in.insert('1.0', 'ATGGCAAGCTTAGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATCGTAGCATCG')
        primer_out = scrolledtext.ScrolledText(t4, bg=Theme.c('bg'), fg=Theme.c('text'),
                                                font=('Courier', 9), height=13, relief='flat', state='disabled')
        primer_out.pack(fill='both', expand=True, padx=12, pady=6)

        def run_primers():
            seq = primer_in.get('1.0', 'end').strip()
            if len(seq) < 40:
                return
            r = design_primers_pro(seq)
            o = primer_out; o.config(state='normal'); o.delete('1.0', 'end')
            o.tag_config('hdr', foreground='#0969da', font=('Courier', 9, 'bold'))
            o.tag_config('good', foreground='#1a7f37')
            o.tag_config('warn', foreground='#9a6700')
            o.tag_config('bad', foreground='#cf222e')
            if 'error' in r:
                o.insert('end', r['error'] + '\n', 'bad')
                o.config(state='disabled'); return
            o.insert('end', ('ПРЯМОЙ праймер (forward):\n' if is_ru else 'FORWARD primer:\n'), 'hdr')
            o.insert('end', '  ' + r['forward']['primer'] + '\n', 'good')
            o.insert('end', '  Tm={}°C  GC={}%  len={}\n\n'.format(
                r['forward']['tm'], r['forward']['gc'], r['forward']['length']))
            o.insert('end', ('ОБРАТНЫЙ праймер (reverse):\n' if is_ru else 'REVERSE primer:\n'), 'hdr')
            o.insert('end', '  ' + r['reverse']['primer'] + '\n', 'good')
            o.insert('end', '  Tm={}°C  GC={}%  len={}\n\n'.format(
                r['reverse']['tm'], r['reverse']['gc'], r['reverse']['length']))
            o.insert('end', ('Контроль качества:\n' if is_ru else 'Quality control:\n'), 'hdr')
            dtm = r['tm_difference']
            o.insert('end', '  ΔTm = {}°C\n'.format(dtm), 'good' if dtm < 3 else 'warn')
            dimer = r['primer_dimer']['risk']
            o.insert('end', '  ' + ('Димер: ' if is_ru else 'Dimer: ') + dimer + '\n',
                      'good' if dimer == 'low' else 'warn' if dimer == 'medium' else 'bad')
            fh = r['forward_hairpin']['risk']
            o.insert('end', '  ' + ('Шпилька F: ' if is_ru else 'Hairpin F: ') + fh + '\n',
                      'good' if fh == 'low' else 'warn' if fh == 'medium' else 'bad')
            o.insert('end', '  ' + ('Размер продукта: ' if is_ru else 'Product size: ')
                      + str(r['product_size']) + ' bp\n')
            o.config(state='disabled')
        tk.Button(t4, text=('▶ Дизайн' if is_ru else '▶ Design'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_primers).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 5: GC Profile + CAI ═══
        t5 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t5, text='GC / CAI')
        tk.Label(t5, text=('GC-профиль, GC-skew (поиск ori), CAI (экспрессия)'
                              if is_ru else 'GC profile, GC-skew (find ori), CAI'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        gc_in = scrolledtext.ScrolledText(t5, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                          font=('Courier', 9), height=4, relief='flat')
        gc_in.pack(fill='x', padx=12, pady=4)
        gc_in.insert('1.0', 'ATGGCGAAACTGGGCGCAGCT' * 8 + 'TTTAAACCCGGG' * 4)
        gc_fig = tk.Frame(t5, bg=Theme.c('bg')); gc_fig.pack(fill='both', expand=True, padx=12, pady=4)
        gc_info = tk.Label(t5, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        gc_info.pack()
        gc_ch = {'c': None}

        def run_gc():
            seq = gc_in.get('1.0', 'end').strip()
            if len(seq) < 50:
                return
            prof = gc_content_profile(seq, window=30)
            skew = gc_skew(seq, window=30)
            cai = codon_adaptation_index(seq)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if gc_ch['c']: gc_ch['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 3.5), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(prof['positions'], prof['gc_profile'], color='#0969da', label='GC%')
            ax.axhline(prof['overall_gc'], color='#57606a', linestyle='--', alpha=0.5)
            ax2 = ax.twinx()
            ax2.plot(skew['positions'], skew['cumulative'], color='#cf222e',
                     label='GC-skew cum', alpha=0.7)
            ax.set_xlabel(('Позиция' if is_ru else 'Position'), color='#1c2128')
            ax.set_ylabel('GC %', color='#0969da')
            ax2.set_ylabel('GC-skew', color='#cf222e')
            ax.tick_params(colors='#57606a'); ax2.tick_params(colors='#57606a')
            ax.grid(alpha=0.2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=gc_fig); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); gc_ch['c'] = cv
            gc_info.config(text='GC={}% | CAI={} ({})'.format(
                prof['overall_gc'], cai['cai'], cai['interpretation']))
        tk.Button(t5, text=('▶ Анализ' if is_ru else '▶ Analyze'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_gc).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 6: Cloning ═══
        t6 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t6, text=('Клонирование' if is_ru else 'Cloning'))
        tk.Label(t6, text=('Gibson assembly — сборка фрагментов по перекрытиям'
                              if is_ru else 'Gibson assembly — join fragments by overlaps'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        clone_in = scrolledtext.ScrolledText(t6, bg=Theme.c('panel2'), fg=Theme.c('text'),
                                              font=('Courier', 9), height=5, relief='flat')
        clone_in.pack(fill='x', padx=12, pady=4)
        clone_in.insert('1.0', 'ATGGCAAGCTTAGCTAGCTACGTGCA\nGCTACGTGCATTTGGGCCCAAAGGGTTT\nGGGTTTACGTACGTAGCTAGCAT')
        clone_out = scrolledtext.ScrolledText(t6, bg=Theme.c('bg'), fg=Theme.c('text'),
                                               font=('Courier', 9), height=11, relief='flat', state='disabled')
        clone_out.pack(fill='both', expand=True, padx=12, pady=6)

        def run_clone():
            frags = [s.strip() for s in clone_in.get('1.0', 'end').split('\n') if s.strip()]
            if len(frags) < 2:
                return
            r = gibson_assembly_simple(frags, overlap=8)
            o = clone_out; o.config(state='normal'); o.delete('1.0', 'end')
            o.tag_config('hdr', foreground='#0969da', font=('Courier', 9, 'bold'))
            o.tag_config('good', foreground='#1a7f37')
            o.tag_config('bad', foreground='#cf222e')
            if 'error' in r:
                o.insert('end', r['error'] + '\n', 'bad')
                o.config(state='disabled'); return
            o.insert('end', ('Результат сборки:\n' if is_ru else 'Assembly result:\n'), 'hdr')
            o.insert('end', '  ' + ('Стыков: ' if is_ru else 'Junctions: ') + str(r['junctions']) + '\n')
            o.insert('end', '  ' + ('Фрагментов: ' if is_ru else 'Fragments: ')
                      + str(r['fragments_used']) + '\n')
            o.insert('end', '  ' + ('Полная сборка: ' if is_ru else 'Complete: '),
                      'hdr')
            o.insert('end', ('да' if r['complete'] else 'нет') if is_ru else
                      ('yes' if r['complete'] else 'no') + '\n',
                      'good' if r['complete'] else 'bad')
            o.insert('end', '\n' + ('Конструкт ({} bp):\n'.format(r['length'])
                                      if is_ru else 'Construct ({} bp):\n'.format(r['length'])), 'hdr')
            assembled = r['assembled']
            for i in range(0, len(assembled), 50):
                o.insert('end', assembled[i:i + 50] + '\n', 'good')
            o.config(state='disabled')
        tk.Button(t6, text=('▶ Собрать' if is_ru else '▶ Assemble'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_clone).pack(fill='x', padx=12, pady=4)

        run_blast()


    def _open_pharmacology_pro(self):
        """Professional pharmacology — PBPK, bioavailability, dose-response."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Фармакология PRO' if is_ru else 'Pharmacology PRO'))
        dlg.geometry('960x720'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Профессиональная фармакология'
                                       if is_ru else 'Professional Pharmacology'),
                  bg=Theme.c('bg'), fg='#8250df', font=('Arial', 15, 'bold')).pack(pady=8)

        nb = ttk.Notebook(dlg); nb.pack(fill='both', expand=True, padx=10, pady=6)

        # ═══ TAB 1: PBPK ═══
        t1 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t1, text='PBPK')
        tk.Label(t1, text=('Физиологическая модель — распределение по органам'
                              if is_ru else 'Physiological model — organ distribution'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        c1 = tk.Frame(t1, bg=Theme.c('bg')); c1.pack(fill='x', padx=12, pady=4)
        p1 = {}
        for i, (k, lbl, d) in enumerate([('dose', 'Доза мг' if is_ru else 'Dose mg', 100),
                ('cl', 'Клиренс L/ч' if is_ru else 'Clearance L/h', 20),
                ('ka', 'ka (всасыв.)' if is_ru else 'ka', 1.0)]):
            r = tk.Frame(c1, bg=Theme.c('bg')); r.pack(side='left', padx=10)
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=7, relief='flat').pack(side='left')
            p1[k] = v
        f1 = tk.Frame(t1, bg=Theme.c('bg')); f1.pack(fill='both', expand=True, padx=12, pady=4)
        info1 = tk.Label(t1, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info1.pack()
        ch1 = {'c': None}

        def run_pbpk():
            dose = safe_float(p1['dose'].get(), 100)
            cl = safe_float(p1['cl'].get(), 20)
            ka = safe_float(p1['ka'].get(), 1.0, minimum=0.01, maximum=5)
            r = pbpk_model(dose, cl, days=1, ka=ka)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch1['c']: ch1['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['blood'], color='#cf222e', label=('Кровь' if is_ru else 'Blood'), linewidth=2)
            colors = {'liver': '#1a7f37', 'kidney': '#0969da', 'brain': '#9a6700',
                      'fat': '#8250df', 'muscle': '#0550ae'}
            for organ, color in colors.items():
                if organ in r['organs']:
                    label = {'liver': 'печень' if is_ru else 'liver',
                             'kidney': 'почки' if is_ru else 'kidney',
                             'brain': 'мозг' if is_ru else 'brain',
                             'fat': 'жир' if is_ru else 'fat',
                             'muscle': 'мышцы' if is_ru else 'muscle'}[organ]
                    ax.plot(r['times'], r['organs'][organ], color=color, label=label, alpha=0.7)
            ax.set_xlabel(('Время (ч)' if is_ru else 'Time (h)'), color='#1c2128')
            ax.set_ylabel(('Концентрация мг/л' if is_ru else 'Concentration mg/L'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=7, ncol=2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f1); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch1['c'] = cv
            info1.config(text='{} {} мг/л {} {} ч'.format(
                ('Пик в крови:' if is_ru else 'Peak blood:'), r['peak_blood_conc'],
                ('на' if is_ru else 'at'), r['peak_time_h']))
        tk.Button(t1, text=('▶ Симулировать' if is_ru else '▶ Simulate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_pbpk).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 2: Bioavailability ═══
        t2 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t2, text=('Биодоступность' if is_ru else 'Bioavailability'))
        tk.Label(t2, text=('Расчёт оральной биодоступности (F = Fa·Fg·Fh)'
                              if is_ru else 'Oral bioavailability (F = Fa·Fg·Fh)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        c2 = tk.Frame(t2, bg=Theme.c('bg')); c2.pack(fill='x', padx=20, pady=8)
        p2 = {}
        for k, lbl, d in [('fa', 'Всасывание Fa (0-1)' if is_ru else 'Absorbed Fa', 0.9),
                          ('fg', 'Кишечник Fg (0-1)' if is_ru else 'Gut Fg', 0.95),
                          ('he', 'Печёночная экстракция (0-1)' if is_ru else 'Hepatic extraction', 0.3)]:
            r = tk.Frame(c2, bg=Theme.c('bg')); r.pack(fill='x', pady=4)
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 10),
                      width=30, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=8, relief='flat').pack(side='left')
            p2[k] = v
        out2 = scrolledtext.ScrolledText(t2, bg=Theme.c('bg'), fg=Theme.c('text'),
                                         font=('Arial', 11), height=10, relief='flat', state='disabled')
        out2.pack(fill='both', expand=True, padx=20, pady=8)

        def run_bioavail():
            fa = safe_float(p2['fa'].get(), 0.9, 0, 1)
            fg = safe_float(p2['fg'].get(), 0.95, 0, 1)
            he = safe_float(p2['he'].get(), 0.3, 0, 1)
            r = bioavailability(fa, fg, he)
            o = out2; o.config(state='normal'); o.delete('1.0', 'end')
            o.tag_config('big', foreground='#1a7f37', font=('Arial', 16, 'bold'))
            o.tag_config('hdr', foreground='#0969da', font=('Arial', 11, 'bold'))
            o.insert('end', ('Биодоступность F:\n' if is_ru else 'Bioavailability F:\n'), 'hdr')
            o.insert('end', '  {}%\n\n'.format(r['bioavailability_pct']), 'big')
            o.insert('end', ('Компоненты:\n' if is_ru else 'Components:\n'), 'hdr')
            o.insert('end', '  Fa (всасывание) = {}\n'.format(r['Fa']))
            o.insert('end', '  Fg (кишечник) = {}\n'.format(r['Fg']))
            o.insert('end', '  Fh (печень) = {}\n'.format(r['Fh']))
            o.insert('end', '\n  ' + ('Потеря при first-pass: ' if is_ru else 'First-pass loss: ')
                      + '{}%\n'.format(r['first_pass_loss_pct']))
            o.config(state='disabled')
        tk.Button(t2, text=('▶ Рассчитать' if is_ru else '▶ Calculate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_bioavail).pack(fill='x', padx=20, pady=4)

        # ═══ TAB 3: Dose-Response ═══
        t3 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t3, text=('Доза-эффект' if is_ru else 'Dose-Response'))
        tk.Label(t3, text=('Кривая доза-эффект (Emax / Hill) + терапевтический индекс'
                              if is_ru else 'Dose-response curve (Emax/Hill) + therapeutic index'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        c3 = tk.Frame(t3, bg=Theme.c('bg')); c3.pack(fill='x', padx=12, pady=4)
        p3d = {}
        for k, lbl, d in [('emax', 'Emax', 100), ('ec50', 'EC50', 10), ('hill', 'Hill', 1.5)]:
            r = tk.Frame(c3, bg=Theme.c('bg')); r.pack(side='left', padx=10)
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=7, relief='flat').pack(side='left')
            p3d[k] = v
        f3 = tk.Frame(t3, bg=Theme.c('bg')); f3.pack(fill='both', expand=True, padx=12, pady=4)
        ch3 = {'c': None}

        def run_dr():
            emax = safe_float(p3d['emax'].get(), 100)
            ec50 = safe_float(p3d['ec50'].get(), 10, minimum=0.01)
            hill = safe_float(p3d['hill'].get(), 1.5, minimum=0.1, maximum=10)
            r = dose_response_curve(emax, ec50, hill)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch3['c']: ch3['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.semilogx(r['doses'], r['effects'], color='#8250df', linewidth=2)
            ax.axhline(emax / 2, color='#57606a', linestyle='--', alpha=0.5)
            ax.axvline(ec50, color='#cf222e', linestyle='--', alpha=0.5, label='EC50')
            ax.set_xlabel(('Концентрация (log)' if is_ru else 'Concentration (log)'), color='#1c2128')
            ax.set_ylabel(('Эффект %' if is_ru else 'Effect %'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f3); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch3['c'] = cv
        tk.Button(t3, text=('▶ Построить' if is_ru else '▶ Plot'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_dr).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 4: Clearance & Half-life ═══
        t4 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t4, text=('Клиренс' if is_ru else 'Clearance'))
        tk.Label(t4, text=('Печёночный клиренс, период полувыведения, терапевт. индекс'
                              if is_ru else 'Hepatic clearance, half-life, therapeutic index'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        c4 = tk.Frame(t4, bg=Theme.c('bg')); c4.pack(fill='x', padx=16, pady=6)
        p4d = {}
        for k, lbl, d in [('clint', 'Внутр. клиренс L/ч' if is_ru else 'Intrinsic CL', 200),
                          ('fu', 'Свободная фракция' if is_ru else 'Free fraction', 0.8),
                          ('vol', 'Объём распред. L' if is_ru else 'Volume L', 50),
                          ('td50', 'TD50 (токсич.)' if is_ru else 'TD50', 500),
                          ('ed50', 'ED50 (эффект.)' if is_ru else 'ED50', 50)]:
            r = tk.Frame(c4, bg=Theme.c('bg')); r.pack(fill='x', pady=3)
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9),
                      width=24, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=8, relief='flat').pack(side='left')
            p4d[k] = v
        out4 = scrolledtext.ScrolledText(t4, bg=Theme.c('bg'), fg=Theme.c('text'),
                                         font=('Arial', 10), height=10, relief='flat', state='disabled')
        out4.pack(fill='both', expand=True, padx=16, pady=8)

        def run_clearance():
            clint = safe_float(p4d['clint'].get(), 200)
            fu = safe_float(p4d['fu'].get(), 0.8, 0, 1)
            vol = safe_float(p4d['vol'].get(), 50, minimum=1)
            td50 = safe_float(p4d['td50'].get(), 500)
            ed50 = safe_float(p4d['ed50'].get(), 50, minimum=0.1)
            hc = hepatic_clearance(clint, 90, fu)
            hl = half_life(hc['hepatic_clearance'], vol)
            ti = therapeutic_index(td50, ed50)
            o = out4; o.config(state='normal'); o.delete('1.0', 'end')
            o.tag_config('hdr', foreground='#0969da', font=('Arial', 10, 'bold'))
            o.tag_config('good', foreground='#1a7f37')
            o.tag_config('bad', foreground='#cf222e')
            o.insert('end', ('Печёночный клиренс:\n' if is_ru else 'Hepatic clearance:\n'), 'hdr')
            o.insert('end', '  CL_h = {} L/ч | E = {}\n'.format(hc['hepatic_clearance'], hc['extraction_ratio']))
            o.insert('end', '  ' + ('тип: ' if is_ru else 'type: ')
                      + ('поток-лимитир.' if hc['flow_limited'] else 'ёмкость-лимитир.' if hc['capacity_limited'] else 'промежуточн.') + '\n\n')
            o.insert('end', ('Период полувыведения:\n' if is_ru else 'Half-life:\n'), 'hdr')
            o.insert('end', '  t½ = {} ч\n'.format(hl['half_life_h']))
            o.insert('end', '  ' + ('до равновесия: ' if is_ru else 'to steady-state: ')
                      + '{} ч\n\n'.format(hl['time_to_steady_state_h']))
            o.insert('end', ('Терапевтический индекс:\n' if is_ru else 'Therapeutic index:\n'), 'hdr')
            o.insert('end', '  TI = {}\n'.format(ti['therapeutic_index']),
                      'good' if not ti['narrow_window'] else 'bad')
            o.insert('end', '  ' + ti['safety'] + '\n',
                      'good' if not ti['narrow_window'] else 'bad')
            o.config(state='disabled')
        tk.Button(t4, text=('▶ Рассчитать' if is_ru else '▶ Calculate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=run_clearance).pack(fill='x', padx=16, pady=4)

        run_pbpk()


    def _open_epidemiology_advanced(self):
        """Advanced epidemiology — SEIRS, data fitting, networks, strains, hospitals."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Эпидемиология PRO+' if is_ru else 'Epidemiology PRO+'))
        dlg.geometry('970x730'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Эпидемиология — продвинутые модели'
                                       if is_ru else 'Epidemiology — advanced models'),
                  bg=Theme.c('bg'), fg='#cf222e', font=('Arial', 15, 'bold')).pack(pady=8)

        nb = ttk.Notebook(dlg); nb.pack(fill='both', expand=True, padx=10, pady=6)

        # ═══ TAB 1: SEIRS (waning immunity) ═══
        t1 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t1, text='SEIRS')
        tk.Label(t1, text=('SEIRS — с потерей иммунитета (повторные волны, как грипп)'
                              if is_ru else 'SEIRS — waning immunity (recurring waves)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        c1 = tk.Frame(t1, bg=Theme.c('bg')); c1.pack(fill='x', padx=12, pady=4)
        p1 = {}
        for k, lbl, d in [('beta', 'β', 0.4), ('immunity', 'Иммунитет (дни)' if is_ru else 'Immunity (days)', 100),
                          ('days', 'Дней' if is_ru else 'Days', 365)]:
            r = tk.Frame(c1, bg=Theme.c('bg')); r.pack(side='left', padx=10)
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=7, relief='flat').pack(side='left')
            p1[k] = v
        f1 = tk.Frame(t1, bg=Theme.c('bg')); f1.pack(fill='both', expand=True, padx=12, pady=4)
        info1 = tk.Label(t1, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info1.pack()
        ch1 = {'c': None}

        def run_seirs():
            beta = safe_float(p1['beta'].get(), 0.4)
            imm = safe_float(p1['immunity'].get(), 100, minimum=1)
            days = safe_int(p1['days'].get(), 365, minimum=10)
            omega = 1.0 / imm
            r = seirs_model(10000, 10, beta, 0.2, 0.1, omega, days)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch1['c']: ch1['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['I'], color='#cf222e', label='I', linewidth=1.5)
            ax.plot(r['times'], r['S'], color='#0969da', label='S', alpha=0.6)
            ax.plot(r['times'], r['R'], color='#1a7f37', label='R', alpha=0.6)
            ax.set_xlabel(('Дни' if is_ru else 'Days'), color='#1c2128')
            ax.set_ylabel(('Люди' if is_ru else 'People'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f1); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch1['c'] = cv
            info1.config(text='R₀={} | {} {} | {} {}'.format(
                r['R0'], ('волн:' if is_ru else 'waves:'), r['waves'],
                ('эндемия' if is_ru else 'endemic') + ':', 'да/yes' if r['endemic'] else 'нет/no'))
        tk.Button(t1, text=('▶ Запустить' if is_ru else '▶ Run'), bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run_seirs).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 2: Data fitting ═══
        t2 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t2, text=('Подгонка' if is_ru else 'Fit data'))
        tk.Label(t2, text=('Найти R₀ из реальных данных заражений'
                              if is_ru else 'Estimate R₀ from real case data'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        tk.Label(t2, text=('Случаи по дням (через запятую):' if is_ru else 'Cases per day (comma-separated):'),
                  bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(anchor='w', padx=12)
        fit_in = tk.Entry(t2, bg=Theme.c('panel2'), fg=Theme.c('text'), relief='flat', font=('Courier', 10))
        fit_in.pack(fill='x', padx=12, pady=4)
        fit_in.insert(0, '10, 16, 25, 40, 63, 100, 158, 251, 398')
        f2 = tk.Frame(t2, bg=Theme.c('bg')); f2.pack(fill='both', expand=True, padx=12, pady=4)
        info2 = tk.Label(t2, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info2.pack()
        ch2 = {'c': None}

        def run_fit():
            try:
                cases = [float(x.strip()) for x in fit_in.get().split(',') if x.strip()]
            except ValueError as _exc:
                biosim_log('warn', 'run_fit: {}'.format(_exc))
                return
            if len(cases) < 3:
                return
            days = list(range(len(cases)))
            r = fit_epidemic_data(days, cases, 100000)
            if 'error' in r:
                info2.config(text=r['error'], fg='#cf222e'); return
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch2['c']: ch2['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.scatter(days, cases, color='#cf222e', label=('Данные' if is_ru else 'Data'), zorder=3)
            ax.plot(r['pred_days'][:len(days) + 10], r['pred_cases'][:len(days) + 10],
                    color='#0969da', label=('Модель' if is_ru else 'Fit'), linestyle='--')
            ax.set_xlabel(('Дни' if is_ru else 'Days'), color='#1c2128')
            ax.set_ylabel(('Случаи' if is_ru else 'Cases'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f2); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch2['c'] = cv
            dt_txt = '{}дн'.format(r['doubling_time_days']) if r['doubling_time_days'] else '—'
            info2.config(text='R₀={} | {} r={} | {} {}'.format(
                r['R0'], ('рост' if is_ru else 'growth'), r['growth_rate'],
                ('удвоение' if is_ru else 'doubling'), dt_txt), fg='#1a7f37')
        tk.Button(t2, text=('▶ Подогнать' if is_ru else '▶ Fit'), bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run_fit).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 3: Network ═══
        t3 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t3, text=('Сеть' if is_ru else 'Network'))
        tk.Label(t3, text=('Распространение по социальной сети (суперспредеры)'
                              if is_ru else 'Spread on social network (superspreaders)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        net_type = tk.StringVar(value='scalefree')
        nrow = tk.Frame(t3, bg=Theme.c('bg')); nrow.pack(pady=4)
        tk.Label(nrow, text=('Тип сети:' if is_ru else 'Network:'), bg=Theme.c('bg'),
                  fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
        for val, lbl in [('scalefree', 'Безмасштабная (с хабами)' if is_ru else 'Scale-free (hubs)'),
                          ('random', 'Случайная' if is_ru else 'Random')]:
            tk.Radiobutton(nrow, text=lbl, variable=net_type, value=val, bg=Theme.c('bg'),
                            fg='#0969da', selectcolor='#eef1f5', font=('Arial', 8)).pack(side='left', padx=4)
        f3 = tk.Frame(t3, bg=Theme.c('bg')); f3.pack(fill='both', expand=True, padx=12, pady=4)
        info3 = tk.Label(t3, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info3.pack()
        ch3 = {'c': None}

        def run_net():
            r = network_epidemic(500, 6, 0.15, 0.1, 100, network_type=net_type.get(), seed=1)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch3['c']: ch3['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(range(len(r['daily_infected'])), r['daily_infected'],
                    color='#cf222e', label=('Заражено' if is_ru else 'Infected'))
            ax.bar(range(len(r['daily_new'])), r['daily_new'], color='#9a6700',
                   alpha=0.4, label=('Новые' if is_ru else 'New'))
            ax.set_xlabel(('Дни' if is_ru else 'Days'), color='#1c2128')
            ax.set_ylabel(('Люди' if is_ru else 'People'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f3); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch3['c'] = cv
            info3.config(text='{} {}% | {} {} | {} {}'.format(
                ('Охват:' if is_ru else 'Attack:'), r['attack_rate'],
                ('пик' if is_ru else 'peak'), r['peak_infected'],
                ('макс.контактов' if is_ru else 'max contacts'), r['max_contacts']))
        tk.Button(t3, text=('▶ Запустить' if is_ru else '▶ Run'), bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run_net).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 4: Strains ═══
        t4 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t4, text=('Штаммы' if is_ru else 'Strains'))
        tk.Label(t4, text=('Конкуренция вариантов вируса (более заразный побеждает)'
                              if is_ru else 'Variant competition (fitter strain wins)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        f4 = tk.Frame(t4, bg=Theme.c('bg')); f4.pack(fill='both', expand=True, padx=12, pady=4)
        info4 = tk.Label(t4, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info4.pack()
        ch4 = {'c': None}

        def run_strains():
            strains = [{'name': 'Альфа' if is_ru else 'Alpha', 'I0': 100, 'beta': 0.3, 'cross_immunity': 0.5},
                       {'name': 'Дельта' if is_ru else 'Delta', 'I0': 10, 'beta': 0.5, 'cross_immunity': 0.5},
                       {'name': 'Омикрон' if is_ru else 'Omicron', 'I0': 1, 'beta': 0.7, 'cross_immunity': 0.3}]
            r = multi_strain_competition(50000, strains, 0.1, 150)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch4['c']: ch4['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            colors = ['#0969da', '#9a6700', '#cf222e']
            for i, name in enumerate(r['names']):
                ax.plot(r['times'], r['strains'][name], color=colors[i % 3], label=name)
            ax.set_xlabel(('Дни' if is_ru else 'Days'), color='#1c2128')
            ax.set_ylabel(('Заражённых' if is_ru else 'Infected'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f4); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch4['c'] = cv
            info4.config(text='{} {} | {}'.format(
                ('Победитель:' if is_ru else 'Winner:'), r['winner'],
                ' / '.join('{}: {}'.format(k, v) for k, v in r['final_sizes'].items())))
        tk.Button(t4, text=('▶ Запустить' if is_ru else '▶ Run'), bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run_strains).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 5: Hospital burden ═══
        t5 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t5, text=('Больницы' if is_ru else 'Hospitals'))
        tk.Label(t5, text=('Нагрузка на больницы и ОРИТ — планирование здравоохранения'
                              if is_ru else 'Hospital/ICU burden — healthcare planning'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        c5 = tk.Frame(t5, bg=Theme.c('bg')); c5.pack(fill='x', padx=12, pady=4)
        p5d = {}
        for k, lbl, d in [('beds', 'Коек' if is_ru else 'Beds', 2000),
                          ('icu', 'ОРИТ' if is_ru else 'ICU', 200)]:
            r = tk.Frame(c5, bg=Theme.c('bg')); r.pack(side='left', padx=10)
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'), width=7, relief='flat').pack(side='left')
            p5d[k] = v
        f5 = tk.Frame(t5, bg=Theme.c('bg')); f5.pack(fill='both', expand=True, padx=12, pady=4)
        info5 = tk.Label(t5, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info5.pack()
        ch5 = {'c': None}

        def run_hospital():
            beds = safe_int(p5d['beds'].get(), 2000)
            icu_cap = safe_int(p5d['icu'].get(), 200)
            sir = sir_model(100000, 100, 0, 0.4, 0.1, 120)
            r = hospital_burden(sir['I'], sir['times'], hospital_capacity=beds, icu_capacity=icu_cap)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch5['c']: ch5['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(sir['times'], r['hospitalized'], color='#9a6700', label=('Госпитализации' if is_ru else 'Hospitalized'))
            ax.plot(sir['times'], r['icu'], color='#cf222e', label='ОРИТ/ICU')
            ax.axhline(beds, color='#9a6700', linestyle='--', alpha=0.5, label=('Коек' if is_ru else 'Bed capacity'))
            ax.axhline(icu_cap, color='#cf222e', linestyle=':', alpha=0.5, label=('ОРИТ макс' if is_ru else 'ICU capacity'))
            ax.set_xlabel(('Дни' if is_ru else 'Days'), color='#1c2128')
            ax.set_ylabel(('Пациентов' if is_ru else 'Patients'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=7)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f5); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch5['c'] = cv
            overflow = ('ПЕРЕПОЛНЕНИЕ' if is_ru else 'OVERFLOW') if r['hospital_overflow'] else ('ОК')
            info5.config(text='{} {} | ОРИТ {} | {}'.format(
                ('Пик коек:' if is_ru else 'Peak beds:'), r['peak_hospitalized'],
                r['peak_icu'], overflow),
                fg='#cf222e' if r['hospital_overflow'] else '#1a7f37')
        tk.Button(t5, text=('▶ Рассчитать' if is_ru else '▶ Calculate'), bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=run_hospital).pack(fill='x', padx=12, pady=4)

        run_seirs()

        # ── organ-system propagation ────────────────────────────────
        # simulate_disease_propagation and get_organ_interactions were both
        # unreachable, so the program held a map of how organ systems drag each
        # other down and never showed it.
        t6e = tk.Frame(nb, bg=Theme.c('bg'))
        nb.add(t6e, text=('Каскад по органам' if is_ru else 'Organ cascade'))
        tk.Label(t6e, text=('Отказ одной системы тянет за собой другие. Здесь '
                            'видно, куда и как быстро расходится нарушение — '
                            'например, почему сердечная недостаточность бьёт по '
                            'почкам, а те по крови.' if is_ru else
                            'Failure in one system drags others down. This shows '
                            'where the disturbance spreads and how fast.'),
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9),
                 wraplength=760, justify='left').pack(anchor='w', padx=12, pady=(6, 0))

        try:
            organ_systems, organ_links = get_organ_interactions()
        except Exception as exc:
            biosim_log('error', 'get_organ_interactions: {}'.format(exc))
            organ_systems, organ_links = {}, []

        ctrl6e = tk.Frame(t6e, bg=Theme.c('bg'))
        ctrl6e.pack(fill='x', padx=12, pady=6)
        tk.Label(ctrl6e, text=('Что отказало:' if is_ru else 'What failed:'),
                 bg=Theme.c('bg'), fg=Theme.c('text'),
                 font=('Arial', 9, 'bold')).pack(side='left')
        sys_keys = list(organ_systems.keys()) if isinstance(organ_systems, dict) \
            else list(organ_systems)

        def sys_label(k):
            v = organ_systems.get(k) if isinstance(organ_systems, dict) else None
            if isinstance(v, dict):
                return v.get('name_ru' if is_ru else 'name_en', k)
            return str(k)

        start_cb = ttk.Combobox(ctrl6e, values=[sys_label(k) for k in sys_keys],
                                state='readonly', font=('Arial', 9), width=28)
        if sys_keys:
            start_cb.current(0)
        start_cb.pack(side='left', padx=6)
        tk.Label(ctrl6e, text=('Шагов:' if is_ru else 'Steps:'),
                 bg=Theme.c('bg'), fg=Theme.c('text'),
                 font=('Arial', 9)).pack(side='left', padx=(10, 2))
        steps_var = tk.StringVar(value='3')
        tk.Entry(ctrl6e, textvariable=steps_var, bg=Theme.c('panel2'),
                 fg=Theme.c('text'), width=5, relief='flat').pack(side='left')

        f6e = tk.Frame(t6e, bg=Theme.c('bg'))
        f6e.pack(fill='both', expand=True, padx=12, pady=6)
        ch6e = {'c': None}
        out6e = scrolledtext.ScrolledText(t6e, bg=Theme.c('panel'),
                                          fg=Theme.c('text'),
                                          font=('Arial', 10), relief='flat',
                                          height=7)
        out6e.pack(fill='x', padx=12, pady=(0, 8))

        def run_cascade():
            if not sys_keys:
                return
            idx = max(0, start_cb.current())
            start_key = sys_keys[idx]
            try:
                n = max(1, min(8, int(safe_float(steps_var.get(), 3))))
            except Exception as exc:
                biosim_log('warn', 'cascade steps: {}'.format(exc))
                n = 3
            try:
                res = simulate_disease_propagation(start_key, n_steps=n)
            except Exception as exc:
                biosim_log('error', 'simulate_disease_propagation: {}'.format(exc))
                out6e.config(state='normal'); out6e.delete('1.0', 'end')
                out6e.insert('end', str(exc)); out6e.config(state='disabled')
                return

            history = res.get('history') if isinstance(res, dict) else res
            if not isinstance(history, list):
                history = [res] if isinstance(res, dict) else []

            out6e.config(state='normal'); out6e.delete('1.0', 'end')
            out6e.tag_config('h', foreground='#0969da', font=('Arial', 10, 'bold'))
            out6e.tag_config('bad', foreground='#cf222e')
            out6e.tag_config('dim', foreground='#57606a', font=('Arial', 9))
            out6e.insert('end', ('Начало: {}\n'.format(sys_label(start_key))
                                 if is_ru else
                                 'Start: {}\n'.format(sys_label(start_key))), 'h')
            final = history[-1] if history else {}
            if isinstance(final, dict):
                ordered = sorted(((k, v) for k, v in final.items()
                                  if isinstance(v, (int, float))),
                                 key=lambda kv: -kv[1])
                for k, v in ordered:
                    tag = 'bad' if v >= 0.5 else None
                    out6e.insert('end', '  {:<28} {:.0%}\n'.format(
                        sys_label(k)[:28], v), tag)
                if len(ordered) > 1:
                    out6e.insert('end', ('\nЗатронуто систем: {} из {}. '
                                         'Чем длиннее каскад, тем труднее лечить '
                                         'первопричину отдельно.\n'
                                         .format(len(ordered), len(sys_keys))
                                         if is_ru else
                                         '\nSystems affected: {} of {}.\n'
                                         .format(len(ordered), len(sys_keys))), 'dim')
            out6e.config(state='disabled')

            try:
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                if ch6e['c']: ch6e['c'].get_tk_widget().destroy()
                fig = plt.Figure(figsize=(8, 3.4), facecolor='#f5f7fa')
                ax = fig.add_subplot(111, facecolor='#ffffff')
                names = set()
                for snap in history:
                    if isinstance(snap, dict):
                        names.update(k for k, v in snap.items()
                                     if isinstance(v, (int, float)))
                colours = ['#cf222e', '#0969da', '#1a7f37', '#9a6700', '#8250df',
                           '#bc4c00', '#0550ae', '#57606a']
                for i, nm in enumerate(sorted(names)):
                    series = []
                    for snap in history:
                        series.append(snap.get(nm, 0.0) * 100
                                      if isinstance(snap, dict) else 0.0)
                    ax.plot(range(len(series)), series,
                            color=colours[i % len(colours)], lw=2,
                            marker='o', markersize=3, label=sys_label(nm)[:20])
                ax.set_xlabel(('Шаг каскада' if is_ru else 'Cascade step'),
                              color='#1c2128')
                ax.set_ylabel(('Нарушение, %' if is_ru else 'Dysfunction, %'),
                              color='#1c2128')
                ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.2)
                ax.tick_params(colors='#57606a')
                for s in ax.spines.values(): s.set_color('#d0d7de')
                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=f6e); cv.draw()
                cv.get_tk_widget().pack(fill='both', expand=True); ch6e['c'] = cv
            except Exception as exc:
                biosim_log('warn', 'run_cascade chart: {}'.format(exc))

        tk.Button(t6e, text=('▶ Проследить каскад' if is_ru else '▶ Trace cascade'),
                  bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                  relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                  pady=4, command=run_cascade).pack(fill='x', padx=12, pady=4)



    def _open_neuroscience_advanced(self):
        """Advanced neuroscience — Izhikevich, cable, Hopfield, perceptron."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Нейробиология PRO+' if is_ru else 'Neuroscience PRO+'))
        dlg.geometry('970x730'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Нейробиология — продвинутые модели'
                                       if is_ru else 'Neuroscience — advanced models'),
                  bg=Theme.c('bg'), fg='#8250df', font=('Arial', 15, 'bold')).pack(pady=8)

        nb = ttk.Notebook(dlg); nb.pack(fill='both', expand=True, padx=10, pady=6)

        # ═══ TAB 1: Izhikevich neuron types ═══
        t1 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t1, text=('Типы нейронов' if is_ru else 'Neuron types'))
        tk.Label(t1, text=('Модель Ижикевича — разные типы реальных нейронов'
                              if is_ru else 'Izhikevich model — real neuron types'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        type_var = tk.StringVar(value='regular_spiking')
        trow = tk.Frame(t1, bg=Theme.c('bg')); trow.pack(pady=4)
        type_opts = [('regular_spiking', 'Регулярный' if is_ru else 'Regular'),
                     ('intrinsically_bursting', 'Пачечный' if is_ru else 'Bursting'),
                     ('chattering', 'Болтающий' if is_ru else 'Chattering'),
                     ('fast_spiking', 'Быстрый' if is_ru else 'Fast')]
        for val, lbl in type_opts:
            tk.Radiobutton(trow, text=lbl, variable=type_var, value=val, bg=Theme.c('bg'),
                            fg='#0969da', selectcolor='#eef1f5', font=('Arial', 8)).pack(side='left', padx=3)
        f1 = tk.Frame(t1, bg=Theme.c('bg')); f1.pack(fill='both', expand=True, padx=12, pady=4)
        info1 = tk.Label(t1, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info1.pack()
        ch1 = {'c': None}

        def run_izh():
            r = izhikevich_neuron(type_var.get(), I_input=10, t_max=200)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch1['c']: ch1['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['V'], color='#8250df', linewidth=0.8)
            ax.set_xlabel(('Время (мс)' if is_ru else 'Time (ms)'), color='#1c2128')
            ax.set_ylabel(('Потенциал (мВ)' if is_ru else 'Potential (mV)'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f1); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch1['c'] = cv
            info1.config(text='{} {} | {}Hz | {}'.format(
                ('Спайков:' if is_ru else 'Spikes:'), r['spike_count'], r['firing_rate_hz'],
                ('пачки' if is_ru else 'bursting') + (': да' if r['bursting'] else ': нет')))
        tk.Button(t1, text=('▶ Симулировать' if is_ru else '▶ Simulate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_izh).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 2: Cable equation ═══
        t2 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t2, text=('Кабель' if is_ru else 'Cable'))
        tk.Label(t2, text=('Распространение сигнала по дендриту (затухание)'
                              if is_ru else 'Signal propagation along dendrite (attenuation)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        f2 = tk.Frame(t2, bg=Theme.c('bg')); f2.pack(fill='both', expand=True, padx=12, pady=4)
        info2 = tk.Label(t2, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info2.pack()
        ch2 = {'c': None}

        def run_cable():
            r = cable_equation(length=50, n_segments=40, I_inject=20, t_max=50)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch2['c']: ch2['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            positions = [i * r['length'] / r['segments'] for i in range(r['segments'])]
            ax.plot(positions, r['final_voltage'], color='#8250df', linewidth=2)
            ax.fill_between(positions, r['final_voltage'], alpha=0.2, color='#8250df')
            ax.set_xlabel(('Расстояние (мкм)' if is_ru else 'Distance (μm)'), color='#1c2128')
            ax.set_ylabel(('Напряжение' if is_ru else 'Voltage'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f2); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch2['c'] = cv
            info2.config(text='{} λ={} | {} τ={}мс'.format(
                ('Постоянная длины' if is_ru else 'Space const'), r['space_constant'],
                ('времени' if is_ru else 'time const'), r['time_constant']))
        tk.Button(t2, text=('▶ Рассчитать' if is_ru else '▶ Calculate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_cable).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 3: Hopfield memory ═══
        t3 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t3, text=('Память' if is_ru else 'Memory'))
        tk.Label(t3, text=('Сеть Хопфилда — ассоциативная память (восстановление образа)'
                              if is_ru else 'Hopfield network — associative memory'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        f3 = tk.Frame(t3, bg=Theme.c('bg')); f3.pack(fill='both', expand=True, padx=12, pady=4)
        info3 = tk.Label(t3, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info3.pack()
        ch3 = {'c': None}

        def run_hopfield():
            # Store 2 patterns (5x5 images), recall noisy version
            # Pattern "T" and "L"
            T = [1,1,1,1,1, -1,-1,1,-1,-1, -1,-1,1,-1,-1, -1,-1,1,-1,-1, -1,-1,1,-1,-1]
            L = [1,-1,-1,-1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1,-1, 1,1,1,1,1]
            # Noisy T (flip 4 pixels)
            noisy = T[:]
            for idx in [2, 7, 12, 20]:
                noisy[idx] *= -1
            r = hopfield_network([T, L], noisy)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch3['c']: ch3['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 3), facecolor='#f5f7fa')
            for i, (title, pat) in enumerate([(('Хранится' if is_ru else 'Stored'), T),
                                               (('Зашумлён' if is_ru else 'Noisy'), noisy),
                                               (('Восстановлен' if is_ru else 'Recalled'), r['recalled'])]):
                ax = fig.add_subplot(1, 3, i + 1, facecolor='#ffffff')
                grid = [[1 if pat[row * 5 + col] > 0 else 0 for col in range(5)] for row in range(5)]
                ax.imshow(grid, cmap='Purples', vmin=0, vmax=1)
                ax.set_title(title, color='#1c2128', fontsize=9)
                ax.set_xticks([]); ax.set_yticks([])
            fig.patch.set_facecolor('#f5f7fa')
            cv = FigureCanvasTkAgg(fig, master=f3); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch3['c'] = cv
            info3.config(text='{} {} {} | {} {}%'.format(
                ('Восстановлено за' if is_ru else 'Converged in'), r['converged_at'],
                ('итер' if is_ru else 'iters'),
                ('качество' if is_ru else 'quality'), r['match_quality']))
        tk.Button(t3, text=('▶ Восстановить образ' if is_ru else '▶ Recall pattern'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_hopfield).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 4: Perceptron learning ═══
        t4 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t4, text=('Обучение' if is_ru else 'Learning'))
        tk.Label(t4, text=('Перцептрон — простейшая обучающаяся нейросеть'
                              if is_ru else 'Perceptron — simplest learning neural net'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        gate_var = tk.StringVar(value='AND')
        grow = tk.Frame(t4, bg=Theme.c('bg')); grow.pack(pady=4)
        tk.Label(grow, text=('Обучить логике:' if is_ru else 'Learn logic:'), bg=Theme.c('bg'),
                  fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
        for g in ['AND', 'OR', 'NAND']:
            tk.Radiobutton(grow, text=g, variable=gate_var, value=g, bg=Theme.c('bg'),
                            fg='#0969da', selectcolor='#eef1f5', font=('Arial', 9)).pack(side='left', padx=4)
        f4 = tk.Frame(t4, bg=Theme.c('bg')); f4.pack(fill='both', expand=True, padx=12, pady=4)
        info4 = tk.Label(t4, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info4.pack()
        ch4 = {'c': None}

        def run_perceptron():
            gate = gate_var.get()
            truth = {'AND': [0,0,0,1], 'OR': [0,1,1,1], 'NAND': [1,1,1,0]}
            inputs = [[0,0],[0,1],[1,0],[1,1]]
            training = list(zip(inputs, truth[gate]))
            r = perceptron_learn(training, epochs=30)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch4['c']: ch4['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(range(1, len(r['accuracy_curve']) + 1),
                    [a * 100 for a in r['accuracy_curve']], color='#1a7f37', marker='o')
            ax.set_xlabel(('Эпоха' if is_ru else 'Epoch'), color='#1c2128')
            ax.set_ylabel(('Точность %' if is_ru else 'Accuracy %'), color='#1c2128')
            ax.set_ylim(0, 105)
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f4); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch4['c'] = cv
            info4.config(text='{} {}% | {} {} | {}'.format(
                ('Точность:' if is_ru else 'Accuracy:'), r['final_accuracy'] * 100,
                ('эпох' if is_ru else 'epochs'), r['epochs_used'],
                ('обучился!' if r['converged'] else 'не сошёлся') if is_ru else
                ('learned!' if r['converged'] else 'not converged')))
        tk.Button(t4, text=('▶ Обучить' if is_ru else '▶ Train'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_perceptron).pack(fill='x', padx=12, pady=4)

        run_izh()

        # ── tuning curve ────────────────────────────────────────────
        # tuning_curve and spike_train_entropy were both in the engine with no
        # way to be reached, so the neuroscience section could simulate a
        # membrane but never ask what a neuron actually encodes.
        t5 = tk.Frame(nb, bg=Theme.c('bg'))
        nb.add(t5, text=('Настройка нейрона' if is_ru else 'Tuning curve'))
        tk.Label(t5, text=('Нейрон отвечает сильнее всего на «свой» стимул и '
                           'слабее на соседние. Ширина настройки решает, что '
                           'важнее: точность или охват.' if is_ru else
                           'A neuron fires hardest at its preferred stimulus and '
                           'less at neighbouring ones. The tuning width decides '
                           'between precision and coverage.'),
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9),
                 wraplength=760, justify='left').pack(anchor='w', padx=12, pady=(6, 0))
        p5n = {}
        c5n = tk.Frame(t5, bg=Theme.c('bg')); c5n.pack(fill='x', padx=12, pady=6)
        for i, (k, lbl, d) in enumerate([
                ('pref', 'Предпочт. угол°' if is_ru else 'Preferred angle°', 90),
                ('width', 'Ширина°' if is_ru else 'Width°', 20),
                ('rate', 'Макс. частота Гц' if is_ru else 'Max rate Hz', 50),
                ('n_cells', 'Нейронов в популяции' if is_ru else 'Cells in population', 6)]):
            r = tk.Frame(c5n, bg=Theme.c('bg'))
            r.grid(row=0, column=i, padx=6, pady=3, sticky='w')
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'),
                     font=('Arial', 9), width=18, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'),
                     width=7, relief='flat').pack(side='left')
            p5n[k] = v
        f5n = tk.Frame(t5, bg=Theme.c('bg'))
        f5n.pack(fill='both', expand=True, padx=12, pady=6)
        ch5n = {'c': None}
        info5n = tk.Label(t5, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                          font=('Arial', 10), wraplength=760, justify='left')
        info5n.pack(fill='x', padx=12)

        def run_tuning():
            width = safe_float(p5n['width'].get(), 20)
            rate = safe_float(p5n['rate'].get(), 50)
            n_cells = max(1, int(safe_float(p5n['n_cells'].get(), 6)))
            stimuli = list(range(0, 181, 2))
            try:
                mine = tuning_curve(safe_float(p5n['pref'].get(), 90), stimuli,
                                    width, rate)['rates']
                pop = []
                for c in range(n_cells):
                    pref = 180.0 * c / n_cells
                    pop.append((pref,
                                tuning_curve(pref, stimuli, width, rate)['rates']))
            except Exception as exc:
                biosim_log('error', 'tuning_curve: {}'.format(exc))
                info5n.config(text=str(exc)); return
            try:
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                if ch5n['c']: ch5n['c'].get_tk_widget().destroy()
                fig = plt.Figure(figsize=(8, 3.4), facecolor='#f5f7fa')
                ax = fig.add_subplot(111, facecolor='#ffffff')
                for pref, curve in pop:
                    ax.plot(stimuli, curve, color='#8250df', alpha=0.3, lw=1)
                ax.plot(stimuli, mine, color='#cf222e', lw=2.5,
                        label=('Выбранный нейрон' if is_ru else 'Selected neuron'))
                ax.set_xlabel(('Ориентация стимула, °' if is_ru
                               else 'Stimulus orientation, °'), color='#1c2128')
                ax.set_ylabel(('Частота разрядов, Гц' if is_ru
                               else 'Firing rate, Hz'), color='#1c2128')
                ax.legend(fontsize=8); ax.grid(alpha=0.2)
                ax.tick_params(colors='#57606a')
                for s in ax.spines.values(): s.set_color('#d0d7de')
                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=f5n); cv.draw()
                cv.get_tk_widget().pack(fill='both', expand=True); ch5n['c'] = cv
            except Exception as exc:
                biosim_log('warn', 'run_tuning chart: {}'.format(exc))
            gap = 180.0 / max(1, n_cells)
            if width < gap * 0.5:
                verdict = ('Настройка узкая, а нейронов мало — между ними «слепые» '
                           'углы, часть стимулов почти никто не кодирует.')
            elif width > gap * 2:
                verdict = ('Настройка широкая — покрытие полное, но соседние углы '
                           'различаются плохо.')
            else:
                verdict = 'Ширина и плотность согласованы — покрытие без провалов.'
            info5n.config(text=verdict if is_ru else
                          'Width {:.0f}° vs spacing {:.0f}°'.format(width, gap))

        tk.Button(t5, text=('▶ Построить' if is_ru else '▶ Plot'), bg='#1a7f37',
                  highlightbackground='#1a7f37', fg='white', relief='flat',
                  font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                  command=run_tuning).pack(fill='x', padx=12, pady=4)

        # ── spike train information ─────────────────────────────────
        t6 = tk.Frame(nb, bg=Theme.c('bg'))
        nb.add(t6, text=('Информация в спайках' if is_ru else 'Spike information'))
        tk.Label(t6, text=('Сколько информации несёт последовательность разрядов. '
                           'Идеально ровный ритм несёт мало информации — как и '
                           'полное молчание.' if is_ru else
                           'How much information a spike train carries. A perfectly '
                           'regular rhythm carries little — and so does silence.'),
                 bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9),
                 wraplength=760, justify='left').pack(anchor='w', padx=12, pady=(6, 0))
        p6n = {}
        c6n = tk.Frame(t6, bg=Theme.c('bg')); c6n.pack(fill='x', padx=12, pady=6)
        for i, (k, lbl, d) in enumerate([
                ('rate', 'Частота Гц' if is_ru else 'Rate Hz', 20),
                ('dur', 'Длительность мс' if is_ru else 'Duration ms', 2000),
                ('bin', 'Окно мс' if is_ru else 'Bin ms', 10),
                ('reg', 'Регулярность 0-1' if is_ru else 'Regularity 0-1', 0.0)]):
            r = tk.Frame(c6n, bg=Theme.c('bg'))
            r.grid(row=0, column=i, padx=6, pady=3, sticky='w')
            tk.Label(r, text=lbl + ':', bg=Theme.c('bg'), fg=Theme.c('text'),
                     font=('Arial', 9), width=16, anchor='w').pack(side='left')
            v = tk.StringVar(value=str(d))
            tk.Entry(r, textvariable=v, bg=Theme.c('panel2'), fg=Theme.c('text'),
                     width=7, relief='flat').pack(side='left')
            p6n[k] = v
        f6n = tk.Frame(t6, bg=Theme.c('bg'))
        f6n.pack(fill='both', expand=True, padx=12, pady=6)
        ch6n = {'c': None}
        info6n = tk.Label(t6, text='', bg=Theme.c('bg'), fg=Theme.c('text'),
                          font=('Arial', 10), wraplength=760, justify='left')
        info6n.pack(fill='x', padx=12)

        def run_spikes():
            import random as _r
            rate = safe_float(p6n['rate'].get(), 20)
            dur = safe_float(p6n['dur'].get(), 2000)
            binsz = safe_float(p6n['bin'].get(), 10)
            reg = max(0.0, min(1.0, safe_float(p6n['reg'].get(), 0.0)))
            _r.seed(7)
            spikes = []
            if rate > 0:
                mean_isi = 1000.0 / rate
                t = 0.0
                while t < dur:
                    jitter = _r.expovariate(1.0 / mean_isi) if reg < 1 else mean_isi
                    isi = mean_isi * reg + jitter * (1 - reg)
                    t += max(1.0, isi)
                    if t < dur:
                        spikes.append(t)
            try:
                r = spike_train_entropy(spikes, dur, binsz)
            except Exception as exc:
                biosim_log('error', 'spike_train_entropy: {}'.format(exc))
                info6n.config(text=str(exc)); return
            try:
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                if ch6n['c']: ch6n['c'].get_tk_widget().destroy()
                fig = plt.Figure(figsize=(8, 3.2), facecolor='#f5f7fa')
                ax = fig.add_subplot(211, facecolor='#ffffff')
                ax.eventplot([spikes], colors=['#0969da'], lineoffsets=0,
                             linelengths=0.8)
                ax.set_xlim(0, dur); ax.set_yticks([])
                ax.set_title(('Растр разрядов' if is_ru else 'Spike raster'),
                             color='#1c2128', fontsize=9)
                for s in ax.spines.values(): s.set_color('#d0d7de')
                ax2 = fig.add_subplot(212, facecolor='#ffffff')
                nb_ = max(1, int(dur / binsz))
                counts = [0] * nb_
                for sp in spikes:
                    idx = min(nb_ - 1, int(sp / binsz))
                    counts[idx] += 1
                ax2.bar([i * binsz for i in range(nb_)], counts, width=binsz * 0.9,
                        color='#1a7f37')
                ax2.set_xlabel(('Время, мс' if is_ru else 'Time, ms'), color='#1c2128')
                ax2.set_ylabel(('Спайков в окне' if is_ru else 'Spikes per bin'),
                               color='#1c2128', fontsize=8)
                ax2.tick_params(colors='#57606a', labelsize=8)
                for s in ax2.spines.values(): s.set_color('#d0d7de')
                fig.tight_layout()
                cv = FigureCanvasTkAgg(fig, master=f6n); cv.draw()
                cv.get_tk_widget().pack(fill='both', expand=True); ch6n['c'] = cv
            except Exception as exc:
                biosim_log('warn', 'run_spikes chart: {}'.format(exc))
            parts = []
            for key in ('entropy_bits', 'mean_rate_hz', 'cv_isi', 'regularity',
                        'bins'):
                if key in r:
                    parts.append('{} = {}'.format(key, r[key]))
            info6n.config(text=' · '.join(parts) if parts else str(r))

        tk.Button(t6, text=('▶ Смоделировать' if is_ru else '▶ Simulate'),
                  bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                  relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                  pady=4, command=run_spikes).pack(fill='x', padx=12, pady=4)



    def _open_synbio_advanced(self):
        """Advanced SynBio — parts library, toggle switch, repressilator, optimization."""
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('SynBio PRO+' if is_ru else 'SynBio PRO+'))
        dlg.geometry('970x730'); dlg.configure(bg=Theme.c('bg'))
        tk.Label(dlg, text=('Синтетическая биология — продвинутые цепи'
                                       if is_ru else 'Synthetic Biology — advanced circuits'),
                  bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 15, 'bold')).pack(pady=8)

        nb = ttk.Notebook(dlg); nb.pack(fill='both', expand=True, padx=10, pady=6)

        # ═══ TAB 1: Parts library ═══
        t1 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t1, text=('Детали' if is_ru else 'Parts'))
        tk.Label(t1, text=('Библиотека биологических деталей (BioBricks/iGEM)'
                              if is_ru else 'Biological parts library (BioBricks/iGEM)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        prow = tk.Frame(t1, bg=Theme.c('bg')); prow.pack(fill='x', padx=12, pady=4)
        tk.Label(prow, text=('Промотор:' if is_ru else 'Promoter:'), bg=Theme.c('bg'),
                  fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
        prom_var = tk.StringVar(value='J23100')
        prom_menu = ttk.Combobox(prow, textvariable=prom_var, state='readonly', width=10,
                                  values=list(list_parts('promoters').keys()))
        prom_menu.pack(side='left', padx=4)
        tk.Label(prow, text='RBS:', bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(side='left', padx=(10, 0))
        rbs_var = tk.StringVar(value='B0034')
        rbs_menu = ttk.Combobox(prow, textvariable=rbs_var, state='readonly', width=8,
                                 values=list(list_parts('rbs').keys()))
        rbs_menu.pack(side='left', padx=4)
        tk.Label(prow, text=('Копий:' if is_ru else 'Copies:'), bg=Theme.c('bg'),
                  fg=Theme.c('text'), font=('Arial', 9)).pack(side='left', padx=(10, 0))
        copy_var = tk.StringVar(value='1')
        tk.Entry(prow, textvariable=copy_var, bg=Theme.c('panel2'), fg=Theme.c('text'), width=5, relief='flat').pack(side='left', padx=4)
        parts_out = scrolledtext.ScrolledText(t1, bg=Theme.c('bg'), fg=Theme.c('text'),
                                               font=('Courier', 9), height=16, relief='flat', state='disabled')
        parts_out.pack(fill='both', expand=True, padx=12, pady=6)

        def run_parts():
            copies = safe_int(copy_var.get(), 1, minimum=1, maximum=100)
            e = calculate_expression(prom_var.get(), rbs_var.get(), copies)
            o = parts_out; o.config(state='normal'); o.delete('1.0', 'end')
            o.tag_config('hdr', foreground='#0969da', font=('Courier', 9, 'bold'))
            o.tag_config('big', foreground='#1a7f37', font=('Courier', 14, 'bold'))
            o.insert('end', ('Расчёт экспрессии:\n' if is_ru else 'Expression calculation:\n'), 'hdr')
            o.insert('end', '  {}%\n\n'.format(e['relative_pct']), 'big')
            o.insert('end', '  ' + ('Промотор ' if is_ru else 'Promoter ') + '{}: {} (сила {})\n'.format(
                e['promoter'], list_parts('promoters')[e['promoter']]['desc'], e['promoter_strength']))
            o.insert('end', '  RBS {}: {} (сила {})\n'.format(
                e['rbs'], list_parts('rbs')[e['rbs']]['desc'], e['rbs_strength']))
            o.insert('end', '  ' + ('Копий: ' if is_ru else 'Copies: ') + str(e['copy_number']) + '\n\n')
            o.insert('end', ('Доступные промоторы:\n' if is_ru else 'Available promoters:\n'), 'hdr')
            for name, p in list_parts('promoters').items():
                o.insert('end', '  {:8s} {:.2f}  {}\n'.format(name, p['strength'], p['desc']))
            o.config(state='disabled')
        tk.Button(t1, text=('▶ Рассчитать' if is_ru else '▶ Calculate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_parts).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 2: Toggle switch ═══
        t2 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t2, text=('Триггер' if is_ru else 'Toggle'))
        tk.Label(t2, text=('Генетический триггер — бистабильная память (2 состояния)'
                              if is_ru else 'Genetic toggle switch — bistable memory'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        srow = tk.Frame(t2, bg=Theme.c('bg')); srow.pack(pady=4)
        tk.Label(srow, text=('Начальное состояние:' if is_ru else 'Initial state:'),
                  bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 9)).pack(side='left')
        init_var = tk.StringVar(value='U')
        tk.Radiobutton(srow, text='U > V', variable=init_var, value='U', bg=Theme.c('bg'),
                        fg='#0969da', selectcolor='#eef1f5', font=('Arial', 9)).pack(side='left', padx=4)
        tk.Radiobutton(srow, text='V > U', variable=init_var, value='V', bg=Theme.c('bg'),
                        fg='#cf222e', selectcolor='#eef1f5', font=('Arial', 9)).pack(side='left', padx=4)
        f2 = tk.Frame(t2, bg=Theme.c('bg')); f2.pack(fill='both', expand=True, padx=12, pady=4)
        info2 = tk.Label(t2, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info2.pack()
        ch2 = {'c': None}

        def run_toggle():
            if init_var.get() == 'U':
                r = toggle_switch(init_u=2.0, init_v=0.5)
            else:
                r = toggle_switch(init_u=0.5, init_v=2.0)
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch2['c']: ch2['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['u'], color='#0969da', label='U (' + ('ген 1' if is_ru else 'gene 1') + ')', linewidth=2)
            ax.plot(r['times'], r['v'], color='#cf222e', label='V (' + ('ген 2' if is_ru else 'gene 2') + ')', linewidth=2)
            ax.set_xlabel(('Время' if is_ru else 'Time'), color='#1c2128')
            ax.set_ylabel(('Концентрация белка' if is_ru else 'Protein concentration'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f2); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch2['c'] = cv
            info2.config(text='{} {} | U={} V={}'.format(
                ('Финальное состояние:' if is_ru else 'Final state:'), r['state'],
                r['final_u'], r['final_v']))
        tk.Button(t2, text=('▶ Симулировать' if is_ru else '▶ Simulate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_toggle).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 3: Repressilator ═══
        t3 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t3, text=('Осциллятор' if is_ru else 'Oscillator'))
        tk.Label(t3, text=('Репрессилятор — генетические часы (3 гена по кругу)'
                              if is_ru else 'Repressilator — genetic clock (3 genes in a ring)'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        f3 = tk.Frame(t3, bg=Theme.c('bg')); f3.pack(fill='both', expand=True, padx=12, pady=4)
        info3 = tk.Label(t3, text='', bg=Theme.c('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        info3.pack()
        ch3 = {'c': None}

        def run_repress():
            r = repressilator()
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if ch3['c']: ch3['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(8, 4), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(r['times'], r['protein1'], color='#0969da', label='LacI')
            ax.plot(r['times'], r['protein2'], color='#cf222e', label='TetR')
            ax.plot(r['times'], r['protein3'], color='#1a7f37', label='cI')
            ax.set_xlabel(('Время' if is_ru else 'Time'), color='#1c2128')
            ax.set_ylabel(('Концентрация белка' if is_ru else 'Protein concentration'), color='#1c2128')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            ax.legend(facecolor='#eef1f5', labelcolor='#1c2128', fontsize=8)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            cv = FigureCanvasTkAgg(fig, master=f3); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); ch3['c'] = cv
            period_txt = '{}'.format(r['period']) if r['period'] else '—'
            info3.config(text='{} {} | {} {} | {} {}'.format(
                ('Колебаний:' if is_ru else 'Oscillations:'), r['oscillations'],
                ('осциллирует' if is_ru else 'oscillating'), 'да/yes' if r['oscillating'] else 'нет/no',
                ('период' if is_ru else 'period'), period_txt))
        tk.Button(t3, text=('▶ Запустить' if is_ru else '▶ Run'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_repress).pack(fill='x', padx=12, pady=4)

        # ═══ TAB 4: Circuit optimization ═══
        t4 = tk.Frame(nb, bg=Theme.c('bg')); nb.add(t4, text=('Оптимизация' if is_ru else 'Optimize'))
        tk.Label(t4, text=('Подбор деталей под нужный уровень экспрессии'
                              if is_ru else 'Find parts for target expression level'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9)).pack(pady=4)
        orow = tk.Frame(t4, bg=Theme.c('bg')); orow.pack(pady=6)
        tk.Label(orow, text=('Целевая экспрессия (0-1):' if is_ru else 'Target expression (0-1):'),
                  bg=Theme.c('bg'), fg=Theme.c('text'), font=('Arial', 10)).pack(side='left')
        target_var = tk.StringVar(value='0.5')
        tk.Entry(orow, textvariable=target_var, bg=Theme.c('panel2'), fg=Theme.c('text'), width=8, relief='flat').pack(side='left', padx=6)
        opt_out = scrolledtext.ScrolledText(t4, bg=Theme.c('bg'), fg=Theme.c('text'),
                                            font=('Courier', 10), height=14, relief='flat', state='disabled')
        opt_out.pack(fill='both', expand=True, padx=12, pady=6)

        def run_optimize():
            target = safe_float(target_var.get(), 0.5, 0, 5)
            r = optimize_circuit(target)
            o = opt_out; o.config(state='normal'); o.delete('1.0', 'end')
            o.tag_config('hdr', foreground='#0969da', font=('Courier', 10, 'bold'))
            o.tag_config('best', foreground='#1a7f37', font=('Courier', 11, 'bold'))
            o.insert('end', ('Цель: ' if is_ru else 'Target: ') + str(target) + '\n\n', 'hdr')
            o.insert('end', ('Лучшая комбинация:\n' if is_ru else 'Best combination:\n'), 'hdr')
            b = r['best_combination']
            o.insert('end', '  {} + {} = {}\n\n'.format(b['promoter'], b['rbs'], b['expression']), 'best')
            o.insert('end', ('Топ-3 варианта:\n' if is_ru else 'Top 3 options:\n'), 'hdr')
            for i, opt in enumerate(r['top_3'], 1):
                o.insert('end', '  {}. {} + {} = {} (откл. {:.3f})\n'.format(
                    i, opt['promoter'], opt['rbs'], opt['expression'], opt['diff']))
            o.config(state='disabled')
        tk.Button(t4, text=('▶ Оптимизировать' if is_ru else '▶ Optimize'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=run_optimize).pack(fill='x', padx=12, pady=4)

        run_parts()


    def _open_ai_assistant(self):
        """Free-form AI assistant powered by Gemini.

        Knows the VIREON program itself (answers "how do I do X here" and
        teaches the user how to use features) AND can search the internet for
        biology / science questions via Gemini's Google Search grounding.
        """
        is_ru = (self.lang == 'ru')
        dlg = tk.Toplevel(self)
        dlg.title('' + ('ИИ-ассистент' if is_ru else 'AI Assistant'))
        dlg.geometry('780x680')
        dlg.configure(bg=Theme.c('bg'))

        tk.Label(dlg, text=('ИИ-ассистент (Gemini)' if is_ru else 'AI Assistant (Gemini)'),
                  bg=Theme.c('bg'), fg='#8250df', font=('Arial', 15, 'bold')).pack(pady=8)
        tk.Label(dlg, text=('Спросите о программе (как что сделать) или любой вопрос по науке — '
                              'по науке ассистент ищет в интернете'
                              if is_ru else
                              'Ask about the program (how to do things) or any science question — '
                              'science questions are looked up online'),
                  bg=Theme.c('bg'), fg=Theme.c('text_dim'), font=('Arial', 9),
                  wraplength=740).pack()

        # Conversation area
        conv = scrolledtext.ScrolledText(dlg, bg=Theme.c('bg'), fg=Theme.c('text'),
                                          font=('Arial', 11), relief='flat',
                                          wrap='word', state='disabled', padx=12, pady=8)
        conv.pack(fill='both', expand=True, padx=14, pady=8)
        conv.tag_config('user', foreground='#0969da', font=('Arial', 11, 'bold'))
        conv.tag_config('ai', foreground='#1a7f37', font=('Arial', 11))
        conv.tag_config('err', foreground='#cf222e', font=('Arial', 11))
        conv.tag_config('dim', foreground=Theme.c('text_dim'), font=('Arial', 9, 'italic'))

        # Conversation history for context
        history = []

        if not ai_is_available():
            conv.config(state='normal')
            conv.insert('end', '' + ai_error_message('no_key', self.lang) + '\n', 'err')
            conv.config(state='disabled')

        # Input area
        inp_frame = tk.Frame(dlg, bg=Theme.c('bg')); inp_frame.pack(fill='x', padx=14, pady=(0, 10))
        question_var = tk.StringVar()
        entry = tk.Entry(inp_frame, textvariable=question_var, bg=Theme.c('panel2'),
                          fg=Theme.c('text'), relief='flat', font=('Arial', 11))
        entry.pack(side='left', fill='x', expand=True, padx=(0, 6), ipady=6)
        entry.focus_set()

        def _is_program_question(text):
            """Heuristic: does the question seem to be about the VIREON program
            itself (how to use it) rather than a general science question?
            Program questions are answered from the system prompt (no web);
            science questions enable internet search."""
            t = text.lower()
            program_words = [
                'biosim', 'программ', 'program', 'app', 'приложени', 'кнопк',
                'button', 'меню', 'menu', 'вкладк', 'tab', 'раздел', 'section',
                'окно', 'window', 'конструктор', 'constructor', 'как сделать',
                'how do i', 'how to use', 'как открыть', 'how to open', 'где найти',
                'where is', 'where can i', 'настройк', 'settings', 'симуляц',
                'simulat', 'канвас', 'canvas', 'экспорт', 'export', 'сохранит',
                'save', 'загрузит', 'load', 'crispr', 'криспр', 'метаболизм',
                'metabolism', 'син.био', 'synbio', 'плазмид', 'plasmid',
                'холст', 'холстов', 'organism', 'организм', 'модул', 'module',
                'шасси', 'chassis', 'библиотек', 'library', 'анализ', 'analysis',
                'дерев', 'tree', 'связ', 'link', 'гайд', 'guide', 'фермент',
                'enzyme', 'нокаут', 'knockout', 'ген ', 'gene', 'стохастик',
                'stochastic', 'монте', 'monte carlo', 'деление клетки',
                'cell division', 'как пользоват', 'how to', 'научи', 'teach me',
            ]
            return any(w in t for w in program_words)

        def send_question(event=None):
            q = question_var.get().strip()
            if not q:
                return
            if not ai_is_available():
                conv.config(state='normal')
                conv.insert('end', '\n' + ai_error_message('no_key', self.lang) + '\n', 'err')
                conv.config(state='disabled'); conv.see('end')
                return
            question_var.set('')
            conv.config(state='normal')
            conv.insert('end', '\n' + ('Вы: ' if is_ru else 'You: '), 'user')
            conv.insert('end', q + '\n', 'user')
            conv.insert('end', ('ИИ думает...' if is_ru else 'AI thinking...') + '\n', 'dim')
            conv.config(state='disabled'); conv.see('end'); dlg.update_idletasks()

            # Watchdog: mark whether a reply arrived; if not within 40s, show a
            # timeout note so the window never stays stuck on "thinking...".
            reply_state = {'done': False}

            def watchdog():
                if reply_state['done'] or not dlg.winfo_exists():
                    return
                conv.config(state='normal')
                try:
                    conv.delete('end-2l', 'end-1l')
                except Exception:
                    pass
                conv.insert('end', ('ИИ не ответил за 40 секунд. Проверьте интернет и ключ в Настройках, попробуйте ещё раз.'
                                     if is_ru else
                                     'No reply within 40s. Check your internet and API key in Settings, then try again.') + '\n', 'err')
                conv.config(state='disabled'); conv.see('end')

            dlg.after(40000, watchdog)

            # Build prompt with history context
            context = ''
            for role, msg in history[-6:]:  # last 6 messages for context
                context += '{}: {}\n'.format(role, msg)
            full_prompt = context + ('Пользователь' if is_ru else 'User') + ': ' + q

            # System prompt: teaches the model about the VIREON program itself.
            program_overview_ru = (
                'СПРАВКА О ПРОГРАММЕ VIREON (симулятор человеческого тела для учёных):\n'
                '• Стартовый экран — выбор языка (рус/англ).\n'
                '• Главное окно: поиск болезней (локальная база + PubMed/базы), '
                'вкладки Инфо, График (ODE-симуляция), 3D-тело, Управление, Молекулы.\n'
                '• Верхнее меню:\n'
                '   – Инструменты: Профиль пациента, ИИ-ассистент (это я), Lab values, '
                'Сравнить болезни, Взаимодействия, Отчёт, Журнал.\n'
                '   – Анализ: PK/PD, Фармакология PRO (PBPK), Фармакокинетика PRO, '
                'Калькулятор доз, Survival (Kaplan-Meier), Статистика, Когорта, Импорт CSV, Manhattan (GWAS).\n'
                '   – Наука: Молекулярка PRO (BLAST, плазмиды), Эпидемиология (SIR/SEIR/SEIRS, сети), '
                'Нейронаука (Ходжкин-Хаксли, сети нейронов), Докинг, MSA, ADME, Викторина.\n'
                '   – Настройки: тут вводится бесплатный ключ Gemini API (aistudio.google.com/apikey) и тема.\n'
                '• Окно СинБио (Синтетическая биология): 4 вкладки сверху — '
                'Конструктор, КРИСПР, Метаболизм, Синтетический организм.\n'
                '\n'
                'КОНСТРУКТОР (сборка генных схем):\n'
                '   – Сверху панель инструментов-категорий: Курсор, Элементы, Белки, '
                'Молекулы, Соединения, Регуляция, Реакции, Входы, Выходы, Модули, '
                'Лаб.данные, Компартменты. Кликаешь категорию → выбираешь часть → '
                'кликаешь на холсте чтобы разместить → соединяешь узлы → '
                'справа кнопка «Запустить симуляцию».\n'
                '   – Виды холста: Схема, Линейная (SBOL), Круг (плазмида), 3D.\n'
                '   – НОВОЕ: несколько холстов. Под строкой вида есть ряд '
                '«Холсты: [Холст 1] [Холст 2] [+]». Кнопка [+] добавляет холст '
                '(до 6), клик по вкладке — переключает (инструменты сверху общие, '
                'меняется только холст), двойной клик — переименовать, правый '
                'клик — удалить. Кнопка «🔗 Связи» открывает диалог, где можно '
                'связать ВЫХОД одного холста со ВХОДОМ другого (сигнал передаётся '
                'между схемами). Кнопка «▶ Симуляция всех» гоняет все холсты '
                'вместе, передавая сигналы по связям — так холсты работают как '
                'связанные подсистемы одного организма.\n'
                '\n'
                'КРИСПР (редактирование генов):\n'
                '   – Слева выбор Организма (человек, E. coli, дрожжи и т.д.) — '
                'организм влияет на геном, лучший фермент Cas, эффективность '
                'HDR/NHEJ, доставку. Кнопка «Профиль организма» показывает геном, '
                'CRISPR-параметры и ключевые гены. (Медицинскую часть про болезни '
                'убрали — CRISPR теперь чистый инструмент редактирования.)\n'
                '   – Загрузка гена по имени (из NCBI) или вставкой '
                'последовательности → выбор фермента → «Найти гайды» → выбор '
                'guide RNA из таблицы → тип правки (нокаут/вставка/замена) → '
                '«Применить». Справа продвинутые инструменты: оценка '
                'эффективности, дизайн донора HDR, предсказание indel, '
                'олигонуклеотиды, праймеры ПЦР, редактор оснований, prime editing, '
                'CRISPRi/CRISPRa, multiplex.\n'
                '\n'
                'МЕТАБОЛИЗМ (модель обмена веществ, стиль COPASI):\n'
                '   – Слева выбор организма и дерево «Модель и задачи». Модель: '
                'виды, реакции, параметры. Задачи: стационар, время, FBA, '
                'бифуркации, Монте-Карло, оптимизация и др. Навигация по дереву '
                'стрелками, открытие Enter.\n'
                '\n'
                'СИНТЕТИЧЕСКИЙ ОРГАНИЗМ (проектирование клетки):\n'
                '   – Слева: выбор Шасси (базовая клетка — E. coli, дрожжи и т.д.) '
                'с её бюджетом ресурсов, кнопка «Экспорт клетки» (CSV генов или '
                'Antimony), и Библиотека модулей в виде дерева: категории '
                '(Метаболизм, Иммунитет, Регуляция, Стресс, Сенсорика...) '
                'раскрываются, модули отмечаются галочкой ✓. Навигация стрелками '
                'вверх/вниз, Enter — включить/выключить модуль или раскрыть '
                'категорию. Есть поиск и «+ Свой модуль».\n'
                '   – Справа вкладки результатов: Обзор, Ресурсы, Рост и '
                'продукция, Геном, Потоки (FBA), Советы.\n'
                '   – НОВОЕ, вкладка «🧬 Модель» (модель клетки во времени): '
                'дерево слева — Гены (промотор, RBS, период полураспада мРНК и '
                'белка), Хромосома, Процессы (число рибосом/РНК-полимераз), '
                'Окружение. В генах есть кнопка «📥 Из CRISPR» (взять '
                'отредактированный ген) и «📥 Из Метаболизма» в Процессах (взять '
                'скорость роста организма).\n'
                '   – НОВОЕ, вкладка «📊 Анализ»: дерево слева с 5 симуляциями — '
                'Симуляция во времени (транскрипция→трансляция→деградация с общим '
                'лимитом рибосом/РНКП), Деление клетки (рост массы + репликация + '
                'деление), Стохастика (шум экспрессии, алгоритм Гиллеспи), '
                'Монте-Карло (разброс по популяции клеток, CV), Чувствительность '
                'роста (скан промотора/RBS → скорость деления). У каждой задачи '
                'кнопка «Запустить», график и числа.\n'
                '• Окно МолБио: работа с последовательностями (BLAST, ORF, рестрикция, дизайн праймеров, плазмиды).\n'
            )
            program_overview_en = (
                'ABOUT THE VIREON PROGRAM (human body simulator for scientists):\n'
                '• Start screen — pick language (RU/EN).\n'
                '• Main window: disease search (local DB + PubMed/online), tabs Info, '
                'Graph (ODE simulation), 3D body, Controls, Molecules.\n'
                '• Top menu:\n'
                '   – Tools: Patient Profile, AI Assistant (that is me), Lab values, '
                'Compare, Interactions, Report, Notes.\n'
                '   – Analysis: PK/PD, Pharmacology PRO (PBPK), Advanced PK, Dose Calculator, '
                'Survival (Kaplan-Meier), Statistics, Cohort, Import CSV, Manhattan (GWAS).\n'
                '   – Science: MolBio PRO (BLAST, plasmids), Epidemiology (SIR/SEIR/SEIRS, networks), '
                'Neuroscience (Hodgkin-Huxley, neural nets), Docking, MSA, ADME, Quiz.\n'
                '   – Settings: enter your free Gemini API key (aistudio.google.com/apikey) and theme here.\n'
                '• SynBio window: 4 top tabs — Constructor, CRISPR, Metabolism, '
                'Synthetic organism.\n'
                '\n'
                'CONSTRUCTOR (building genetic circuits):\n'
                '   – Top toolbar of category buttons: Cursor, Elements, Proteins, '
                'Molecules, Connections, Regulation, Reactions, Inputs, Outputs, '
                'Modules, Lab data, Compartments. Click a category → pick a part → '
                'click on the canvas to place it → connect nodes → "Run simulation" '
                'on the right.\n'
                '   – Canvas views: Diagram, Linear (SBOL), Circular (plasmid), 3D.\n'
                '   – NEW: multiple canvases. Below the view row there is a row '
                '"Canvases: [Canvas 1] [Canvas 2] [+]". [+] adds a canvas (up to '
                '6), clicking a tab switches (the top tools are shared, only the '
                'canvas changes), double-click renames, right-click deletes. The '
                '"🔗 Links" button opens a dialog to link an OUTPUT of one canvas '
                "to an INPUT of another (a signal passes between circuits). "
                '"▶ Simulate all" runs every canvas together, passing signals '
                'along the links, so canvases act as coupled subsystems of one '
                'organism.\n'
                '\n'
                'CRISPR (gene editing):\n'
                '   – On the left, pick an Organism (human, E. coli, yeast, etc.) — '
                'the organism sets the genome, best Cas enzyme, HDR/NHEJ '
                'efficiency and delivery. "Organism profile" shows the genome, '
                'CRISPR parameters and key genes. (The medical/disease part was '
                'removed — CRISPR is now a pure editing tool.)\n'
                '   – Load a gene by name (from NCBI) or by pasting a sequence → '
                'pick an enzyme → "Find guides" → choose a guide RNA from the '
                'table → edit type (knockout/insertion/replacement) → "Apply". '
                'Advanced tools on the right: efficiency scoring, HDR donor '
                'design, indel prediction, cloning oligos, PCR primers, base '
                'editor, prime editing, CRISPRi/CRISPRa, multiplex.\n'
                '\n'
                'METABOLISM (metabolic model, COPASI-style):\n'
                '   – Left: organism picker and a "Model & Tasks" tree. Model: '
                'species, reactions, parameters. Tasks: steady state, time '
                'course, FBA, bifurcation, Monte Carlo, optimization, etc. '
                'Navigate the tree with arrows, open with Enter.\n'
                '\n'
                'SYNTHETIC ORGANISM (designing a cell):\n'
                '   – Left: pick a Chassis (base cell — E. coli, yeast, etc.) with '
                'its resource budget, an "Export cell" button (genes CSV or '
                'Antimony), and a Module Library as a tree: categories '
                '(Metabolism, Immunity, Regulation, Stress, Sensing...) expand, '
                'modules are checked with ✓. Navigate with Up/Down arrows, Enter '
                'toggles a module or opens a category. There is search and '
                '"+ Custom module".\n'
                '   – Right, result tabs: Overview, Resources, Growth & output, '
                'Genome, Fluxes (FBA), Advice.\n'
                '   – NEW "🧬 Model" tab (time-resolved cell model): a left tree — '
                'Genes (promoter, RBS, mRNA & protein half-life), Chromosome, '
                'Processes (ribosome/RNA-polymerase counts), Environment. Genes '
                'has a "📥 From CRISPR" button (pull an edited gene) and Processes '
                'has "📥 From Metabolism" (pull the organism growth rate).\n'
                '   – NEW "📊 Analysis" tab: a left tree of 5 simulations — Time '
                'Course (transcription→translation→degradation with a shared '
                'ribosome/RNAP limit), Cell Division (mass growth + replication + '
                'division), Stochastic (expression noise, Gillespie), Monte Carlo '
                '(cell-to-cell spread, CV), Growth Sensitivity (scan promoter/RBS '
                '→ division speed). Each task has a Run button, a chart and '
                'numbers.\n'
                '• MolBio window: sequence work (BLAST, ORF, restriction, primer design, plasmids).\n'
            )
            is_prog = _is_program_question(q)
            if is_ru:
                system = (
                    'Ты — встроенный ассистент образовательной программы VIREON по биологии и медицине. '
                    'Ты знаешь программу и умеешь учить пользователя пользоваться ей. '
                    'Если вопрос о программе (как что-то сделать, где найти функцию) — объясни пошагово, '
                    'опираясь на справку ниже, и научи, куда нажимать. '
                    'Если вопрос по науке (биология, химия, медицина) — отвечай точно и научно, '
                    'при необходимости используй свежие данные из интернета. Приводи примеры. '
                    'Для медицинских тем напоминай про консультацию врача.\n\n'
                    + program_overview_ru)
            else:
                system = (
                    'You are the built-in assistant of the VIREON biology/medicine education app. '
                    'You know the program and teach the user how to use it. '
                    'If the question is about the program (how to do something, where a feature is) — '
                    'explain step by step using the reference below and teach where to click. '
                    'If the question is about science (biology, chemistry, medicine) — answer accurately '
                    'and scientifically, using fresh info from the internet when helpful. Give examples. '
                    'For medical topics, remind about consulting a doctor.\n\n'
                    + program_overview_en)

            # Enable web search only for science questions, not program help.
            # Run the network call in a background thread so the window does
            # NOT freeze while waiting for the reply (Tkinter is single-threaded;
            # a blocking call here would hang the whole UI).
            def worker():
                # Catch EVERYTHING: in a packaged .app/.exe an uncaught error
                # (e.g. SSL) would kill this thread silently and leave
                # "thinking..." on screen forever. We always deliver a result.
                try:
                    r = call_gemini_api(full_prompt, system_prompt=system,
                                         max_tokens=1500, use_search=(not is_prog))
                except Exception as _e:
                    biosim_log('warn', 'ai_assistant worker: {}'.format(_e))
                    r = {'ok': False, 'text': '',
                         'error': ('Ошибка соединения: {}'.format(_e)[:200] if is_ru
                                    else 'Connection error: {}'.format(_e)[:200])}

                def show_result():
                    reply_state['done'] = True
                    if not dlg.winfo_exists():
                        return  # window was closed while waiting
                    conv.config(state='normal')
                    # Remove "thinking..." line
                    conv.delete('end-2l', 'end-1l')
                    if r['ok']:
                        conv.insert('end', ('ИИ: ' if is_ru else 'AI: '), 'ai')
                        conv.insert('end', r['text'] + '\n', 'ai')
                        history.append(('User', q))
                        history.append(('Assistant', r['text']))
                    else:
                        conv.insert('end', '' + ai_error_message(r['error'], self.lang) + '\n', 'err')
                    conv.config(state='disabled'); conv.see('end')

                # Hand the UI update back to the main Tk thread.
                try:
                    dlg.after(0, show_result)
                except Exception:
                    pass

            threading.Thread(target=worker, daemon=True).start()

        entry.bind('<Return>', send_question)
        tk.Button(inp_frame, text=('Отправить' if is_ru else 'Send'),
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1, font=('Arial', 10, 'bold'),
                   cursor='hand2', padx=16, command=send_question).pack(side='left')

        # Quick example questions
        examples_frame = tk.Frame(dlg, bg=Theme.c('bg')); examples_frame.pack(fill='x', padx=14, pady=(0, 8))
        examples = (['Как запустить симуляцию?', 'Что такое CRISPR?', 'Как работает инсулин?']
                     if is_ru else ['How do I run a simulation?', 'What is CRISPR?', 'How does insulin work?'])
        for ex in examples:
            tk.Button(examples_frame, text=ex, bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'),
                       relief='solid', bd=1, font=('Arial', 8), cursor='hand2',
                       command=lambda e=ex: (question_var.set(e), send_question())
                       ).pack(side='left', padx=2)


    def _open_anatomy_quiz(self):
        """Anatomy quiz — guess organ functions / locations."""
        dlg = tk.Toplevel(self)
        dlg.title('' + ('Викторина по анатомии' if self.lang == 'ru'
                              else 'Anatomy Quiz'))
        dlg.geometry('680x540')
        dlg.configure(bg=Theme.c('bg'))
        dlg.transient(self)

        # Quiz questions
        questions_ru = [
            {'q': 'Какой орган производит инсулин?',
              'opts': ['Печень', 'Поджелудочная железа', 'Селезёнка', 'Почки'],
              'ans': 1, 'hint': 'β-клетки в островках Лангерганса'},
            {'q': 'Где образуются эритроциты?',
              'opts': ['Лёгкие', 'Костный мозг', 'Печень', 'Селезёнка'],
              'ans': 1, 'hint': 'У взрослых — в губчатых костях'},
            {'q': 'Сколько камер в сердце человека?',
              'opts': ['2', '3', '4', '5'],
              'ans': 2, 'hint': '2 предсердия + 2 желудочка'},
            {'q': 'Какой гормон отвечает за стресс?',
              'opts': ['Инсулин', 'Кортизол', 'Тироксин', 'Мелатонин'],
              'ans': 1, 'hint': 'Производится надпочечниками'},
            {'q': 'Где находится мозжечок?',
              'opts': ['Лоб', 'Затылок', 'Висок', 'Темя'],
              'ans': 1, 'hint': 'Под затылочными долями'},
            {'q': 'Какой орган фильтрует кровь?',
              'opts': ['Сердце', 'Печень', 'Почки', 'Селезёнка'],
              'ans': 2, 'hint': 'У человека их два, по обе стороны позвоночника'},
            {'q': 'Сколько хромосом у человека?',
              'opts': ['23', '46', '48', '64'],
              'ans': 1, 'hint': '23 пары = 46 всего'},
            {'q': 'Какая часть мозга отвечает за память?',
              'opts': ['Мозжечок', 'Гиппокамп', 'Гипофиз', 'Гипоталамус'],
              'ans': 1, 'hint': 'Часть лимбической системы'},
        ]
        questions_en = [
            {'q': 'Which organ produces insulin?',
              'opts': ['Liver', 'Pancreas', 'Spleen', 'Kidneys'],
              'ans': 1, 'hint': 'β-cells in islets of Langerhans'},
            {'q': 'Where are red blood cells made?',
              'opts': ['Lungs', 'Bone marrow', 'Liver', 'Spleen'],
              'ans': 1, 'hint': 'In adults — spongy bones'},
            {'q': 'How many chambers in the human heart?',
              'opts': ['2', '3', '4', '5'],
              'ans': 2, 'hint': '2 atria + 2 ventricles'},
            {'q': 'Which hormone is the stress hormone?',
              'opts': ['Insulin', 'Cortisol', 'Thyroxine', 'Melatonin'],
              'ans': 1, 'hint': 'Produced by adrenals'},
            {'q': 'Where is the cerebellum?',
              'opts': ['Forehead', 'Back of head', 'Temple', 'Crown'],
              'ans': 1, 'hint': 'Below occipital lobes'},
            {'q': 'Which organ filters blood?',
              'opts': ['Heart', 'Liver', 'Kidneys', 'Spleen'],
              'ans': 2, 'hint': 'Humans have two, on either side of spine'},
            {'q': 'How many chromosomes do humans have?',
              'opts': ['23', '46', '48', '64'],
              'ans': 1, 'hint': '23 pairs = 46 total'},
            {'q': 'Which brain part handles memory?',
              'opts': ['Cerebellum', 'Hippocampus', 'Pituitary', 'Hypothalamus'],
              'ans': 1, 'hint': 'Part of limbic system'},
        ]

        questions = questions_ru if self.lang == 'ru' else questions_en
        import random
        random.shuffle(questions)
        questions = questions[:5]  # take 5

        state = {'idx': 0, 'score': 0, 'answers': []}

        # Header
        tk.Label(dlg, text=('Викторина по анатомии' if self.lang == 'ru'
                                       else 'Anatomy Quiz'),
                  bg=Theme.c('bg'), fg='#9a6700',
                  font=('Arial', 15, 'bold')).pack(pady=10)

        progress_lbl = tk.Label(dlg, text='', bg=Theme.c('bg'), fg=Theme.c('text_dim'),
                                    font=('Arial', 10))
        progress_lbl.pack()

        # Question card
        q_card = tk.Frame(dlg, bg=Theme.c('panel'),
                             highlightthickness=2,
                             highlightbackground=Theme.c('border'))
        q_card.pack(fill='both', expand=True, padx=20, pady=10)

        q_lbl = tk.Label(q_card, text='', bg=Theme.c('panel'), fg=Theme.c('text'),
                            font=('Arial', 13, 'bold'),
                            wraplength=600, justify='left', padx=14, pady=14)
        q_lbl.pack(anchor='w')

        # Options
        opts_frame = tk.Frame(q_card, bg=Theme.c('panel'))
        opts_frame.pack(fill='both', expand=True, padx=14, pady=8)
        opt_buttons = []
        for i in range(4):
            b = tk.Button(opts_frame, text='',
                            bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                            font=('Arial', 11), cursor='hand2',
                            padx=14, pady=8, anchor='w', justify='left',
                            wraplength=550)
            b.pack(fill='x', pady=4)
            opt_buttons.append(b)

        # Feedback area
        feedback_lbl = tk.Label(q_card, text='', bg=Theme.c('panel'),
                                    font=('Arial', 10, 'italic'),
                                    wraplength=600, justify='left')
        feedback_lbl.pack(padx=14, pady=4)

        def show_question():
            if state['idx'] >= len(questions):
                # Show final score
                q_lbl.config(text=('Викторина завершена!' if self.lang == 'ru'
                                             else 'Quiz complete!'))
                for b in opt_buttons:
                    b.pack_forget()
                fb = '' + ('Правильных ответов: ' if self.lang == 'ru' else 'Correct: ')
                fb += '{}/{}'.format(state['score'], len(questions))
                pct = 100 * state['score'] / len(questions)
                if pct >= 80:
                    fb += '\n' + ('Отличный результат!' if self.lang == 'ru'
                                       else 'Excellent!')
                    feedback_lbl.config(text=fb, fg='#1a7f37')
                elif pct >= 60:
                    fb += '\n' + ('Хорошо!' if self.lang == 'ru' else 'Good job!')
                    feedback_lbl.config(text=fb, fg='#9a6700')
                else:
                    fb += '\n' + ('Стоит повторить материал.' if self.lang == 'ru'
                                      else 'Keep studying!')
                    feedback_lbl.config(text=fb, fg='#bc4c00')
                progress_lbl.config(text='')
                return

            q = questions[state['idx']]
            progress_lbl.config(text='{}/{} · Score: {}'.format(
                state['idx'] + 1, len(questions), state['score']))
            q_lbl.config(text='Q{}: {}'.format(state['idx'] + 1, q['q']))
            feedback_lbl.config(text='')

            for i, opt in enumerate(q['opts']):
                opt_buttons[i].config(text=chr(65 + i) + ') ' + opt,
                                          bg=Theme.c('panel2'), fg=Theme.c('text'),
                                          state='normal',
                                          command=lambda idx=i: answer(idx))
                opt_buttons[i].pack(fill='x', pady=4)

        def answer(idx):
            q = questions[state['idx']]
            correct = (idx == q['ans'])
            if correct:
                state['score'] += 1
                opt_buttons[idx].config(bg='#1a7f37', fg='white')
                feedback_lbl.config(text=('Правильно!' if self.lang == 'ru'
                                                     else 'Correct!') + ' ' + q['hint'],
                                        fg='#1a7f37')
            else:
                opt_buttons[idx].config(bg='#cf222e', fg='white')
                opt_buttons[q['ans']].config(bg='#1a7f37', fg='white')
                feedback_lbl.config(text=('Правильный ответ: ' if self.lang == 'ru'
                                                     else 'Correct answer: ')
                                            + q['opts'][q['ans']] + '. ' + q['hint'],
                                        fg='#cf222e')

            for b in opt_buttons:
                b.config(state='disabled')

            state['idx'] += 1
            dlg.after(1800, show_question)

        show_question()

        tk.Button(dlg, text='Close',
                   bg=Theme.c('panel'), highlightbackground=Theme.c('border'), fg=Theme.c('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy
                   ).pack(fill='x', padx=14, pady=8)


    def _reset(self):
        self.cur_dis = None; self.cur_dat = None
        self.drug_var.set(''); self.drug_dose.set(50)
        self.mut_val.set(0); self.days_var.set(100)
        self.dlbl.config(text='100'); self.dose_lbl.config(text='50')
        self.ax.clear(); self._sax()
        self.ax.text(0.5, 0.5, self.L['hint'], transform=self.ax.transAxes,
                     ha='center', va='center', color='#57606a')
        self.cv.draw()
        self.human.draw_default(self.L['hint'])
        self._iw([(self.L['hint'] + '\n', 'norm')])
        self._rw([(self.L['hint'] + '\n', 'nrm')])


# ─────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────
if __name__ == '__main__':
    LangScreen().mainloop()

