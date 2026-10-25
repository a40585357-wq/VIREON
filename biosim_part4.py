"""
VIREON Part 4 — MolBio Window (molecular biology workspace).
Split from biosim_part2.py to keep file sizes manageable and reduce editor lag.
All shared data/functions stay in part2/part3; this file only holds the
self-contained MolBioWindow UI class.
"""
import os, sys, math, threading, json, random
import urllib.request, urllib.parse, ssl
import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np


# ── Logger ───────────────────────────────────────────────────────────
# biosim_log is really defined in biosim_part2. Import order (and packaged
# .app/.exe builds) can leave it undefined in this module when one of the 44
# calls below runs, raising "name 'biosim_log' is not defined". This wrapper
# looks up the real logger at CALL time and forwards to it, else prints.
def biosim_log(level, message):
    """Forward to biosim_part2.biosim_log if available, else print."""
    import sys as _sys
    _p2 = _sys.modules.get('biosim_part2')
    _real = getattr(_p2, 'biosim_log', None) if _p2 is not None else None
    if _real is not None and _real is not biosim_log:
        try:
            return _real(level, message)
        except Exception:
            pass
    try:
        print('[biosim] {}: {}'.format(level, message))
    except Exception:
        pass


try:
    from biosim_part3 import *
except Exception as _exc:
    # biosim_log is defined in biosim_part2, which may not be imported yet when
    # part4 loads; guard these calls to avoid a NameError during import.
    if 'biosim_log' in globals():
        biosim_log('warn', 'module: {}'.format(_exc))
try:
    from PIL import Image, ImageTk
    HAS_PIL = True
except Exception as _exc:
    if 'biosim_log' in globals():
        biosim_log('warn', 'module: {}'.format(_exc))
    HAS_PIL = False

# SPEED: rdkit is heavy (~5-10s on old machines). Lazy-load it only when a
# chemistry feature is actually used, instead of at startup.
HAS_RDKIT = None
Chem = None
AllChem = Draw = None


def _ensure_rdkit_p4():
    """Load rdkit on first use (part4 chemistry features)."""
    global HAS_RDKIT, Chem, AllChem, Draw
    if HAS_RDKIT is None:
        try:
            from rdkit import Chem as _Chem
            from rdkit.Chem import AllChem as _AllChem, Draw as _Draw
            Chem = _Chem; AllChem = _AllChem; Draw = _Draw
            HAS_RDKIT = True
        except Exception as _exc:
            biosim_log('warn', '_ensure_rdkit_p4: {}'.format(_exc))
            HAS_RDKIT = False
    return HAS_RDKIT



# These molecular-biology helpers live in part2. We import them lazily to avoid
# a circular import (part2 imports MolBioWindow from this file).
def _import_part2_helpers():
    import biosim_part2 as _p2
    g = globals()
    for _name in dir(_p2):
        if _name not in g and not _name.startswith('__'):
            g[_name] = getattr(_p2, _name)

try:
    _import_part2_helpers()
except Exception as _exc:
    biosim_log('warn', 'module: {}'.format(_exc))


class MolBioWindow(tk.Toplevel):
    """Molecular Biology Workspace.
    Sequence editing, primer design, restriction analysis,
    alignment, database lookup, structure prediction.
    
    6 tabs:
        Sequence Editor — DNA/RNA/protein editor with analysis
        Primer Designer — PCR primer design with Tm/GC/hairpins
        Restriction & Cloning — sites, digestion, virtual gel
        Alignment — pairwise + BLAST-like search
        Databases — UniProt/NCBI fetch
        Structure — RNA secondary, codon optimization
    """

    TR = {
        'ru': {
            'title': 'Молекулярная биология — VIREON',
            'tab_seq': 'Последовательности',
            'tab_primer': 'Праймеры',
            'tab_restriction': 'Ферменты и клонирование',
            'tab_align': 'Выравнивание',
            'tab_database': 'Базы данных',
            'tab_structure': 'Структура',
            'tab_crispr': 'CRISPR',
            'tab_phylo': 'Филогения',
            'welcome_msg': 'Добро пожаловать в Молекулярную Биологию',
            'subtitle': 'Sequence editing · Primer design · Cloning · BLAST · Structure',
        },
        'en': {
            'title': 'Molecular Biology — VIREON',
            'tab_seq': 'Sequences',
            'tab_primer': 'Primers',
            'tab_restriction': 'Restriction & Cloning',
            'tab_align': 'Alignment',
            'tab_database': 'Databases',
            'tab_structure': 'Structure',
            'tab_crispr': 'CRISPR',
            'tab_phylo': 'Phylogeny',
            'welcome_msg': 'Welcome to Molecular Biology',
            'subtitle': 'Sequence editing · Primer design · Cloning · BLAST · Structure',
        },
    }

    def __init__(self, parent, lang='ru'):
        _import_part2_helpers()  # load molbio functions from part2
        super().__init__(parent)
        self.parent_main = parent
        self.lang = lang if lang in ('ru', 'en') else 'en'
        self.T = self.TR[self.lang]

        self.title(self.T['title'])
        self.geometry('1200x780')
        self.configure(bg=biosim_theme('bg'))
        self.minsize(1000, 650)

        # Shared state across tabs
        self.current_sequence = ''  # last loaded/edited sequence
        self.sequence_type = 'dna'  # 'dna' | 'rna' | 'protein'
        self.sequence_name = 'Untitled'
        self.history = []  # for undo

        self._build_ui()
        biosim_log('info', 'MolBio Workspace opened')

    def _build_ui(self):
        # ── Top bar ──
        top = tk.Frame(self, bg=biosim_theme('panel'), height=50)
        top.pack(fill='x')
        top.pack_propagate(False)

        tk.Label(top, text='',
                  bg=biosim_theme('panel'), fg='#8250df',
                  font=('Arial', 22)).pack(side='left', padx=12, pady=4)
        title_frame = tk.Frame(top, bg=biosim_theme('panel'))
        title_frame.pack(side='left', pady=4)
        tk.Label(title_frame, text='Molecular Biology',
                  bg=biosim_theme('panel'), fg=biosim_theme('text'),
                  font=('Arial', 14, 'bold')).pack(anchor='w')
        tk.Label(title_frame, text=self.T['subtitle'],
                  bg=biosim_theme('panel'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9, 'italic')).pack(anchor='w')

        # Quick-action buttons in top bar
        quick_frame = tk.Frame(top, bg=biosim_theme('panel'))
        quick_frame.pack(side='right', padx=8)
        for icon, txt, color, cmd in [
            ('', 'Paste',
              '#0969da', self._quick_paste),
            ('', 'Open',
              '#1a7f37', self._quick_open_file),
            ('', 'Fetch',
              '#8250df', self._quick_fetch_ncbi),
        ]:
            btn = tk.Button(quick_frame, text=icon + ' ' + txt,
                              bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=color, relief='solid', bd=1,
                              font=('Arial', 9, 'bold'), cursor='hand2',
                              padx=8, pady=3,
                              command=cmd)
            btn.pack(side='left', padx=2)

        # ── Notebook with tabs ──
        style = ttk.Style()
        try:
            style.theme_use('clam')
        except Exception as _exc:
            biosim_log('warn', '_build_ui: {}'.format(_exc))
        style.configure('MB.TNotebook', background='#f5f7fa',
                          borderwidth=0)
        style.configure('MB.TNotebook.Tab',
                          background='#ffffff', foreground=biosim_theme('text'),
                          padding=(14, 6), borderwidth=0,
                          font=('Arial', 10, 'bold'))
        style.map('MB.TNotebook.Tab',
                    background=[('selected', '#eef1f5')],
                    foreground=[('selected', '#8250df')])

        self.notebook = ttk.Notebook(self, style='MB.TNotebook')
        self.notebook.pack(fill='both', expand=True, padx=8, pady=8)

        # Build all tabs into GROUPED sub-notebooks for a clean, compact layout.
        # Each group is a top-level tab; the real tabs live inside as sub-tabs.
        is_ru = (self.lang == 'ru')
        self._mb_groups = {}

        def make_group(group_key, label):
            frame = tk.Frame(self.notebook, bg=biosim_theme('bg'))
            self.notebook.add(frame, text=label)
            sub_nb = ttk.Notebook(frame, style='MB.TNotebook')
            sub_nb.pack(fill='both', expand=True, padx=2, pady=2)
            self._mb_groups[group_key] = sub_nb
            return sub_nb

        # Create the 5 thematic groups
        grp_seq = make_group('seq', '' + ('Последовательность' if is_ru else 'Sequence'))
        grp_align = make_group('align', '' + ('Выравнивание' if is_ru else 'Alignment'))
        grp_struct = make_group('struct', '' + ('Структура' if is_ru else 'Structure'))
        grp_tools = make_group('tools', '' + ('Инструменты' if is_ru else 'Tools'))
        grp_analysis = make_group('analysis', '' + ('Анализ' if is_ru else 'Analysis'))

        # Tell each _build_tab method which sub-notebook to use by temporarily
        # pointing self.notebook at the right group while building.
        _real_notebook = self.notebook

        def build_into(group_nb, *build_methods):
            self.notebook = group_nb
            for m in build_methods:
                m()
            self.notebook = _real_notebook

        # Sequence group: raw sequence work + DNA analysis
        build_into(grp_seq, self._build_tab_sequence, self._build_tab_orf,
                   self._build_tab_rna, self._build_tab_gcanalysis)
        # Alignment group: all sequence comparison
        build_into(grp_align, self._build_tab_alignment, self._build_tab_blast_pro,
                   self._build_tab_msa, self._build_tab_phylo)
        # Structure group: 3D/structure
        build_into(grp_struct, self._build_tab_structure, self._build_tab_protein)
        # Tools group: cloning & editing tools
        build_into(grp_tools, self._build_tab_primer, self._build_tab_restriction,
                   self._build_tab_digest, self._build_tab_crispr)
        # Analysis group: databases, annotation, drug analysis
        build_into(grp_analysis, self._build_tab_databases, self._build_tab_annotation_pro,
                   self._build_tab_docking, self._build_tab_adme)


        # ── Status bar ──
        self.status_bar = tk.Frame(self, bg=biosim_theme('panel'), height=22)
        self.status_bar.pack(side='bottom', fill='x')
        self.status_bar.pack_propagate(False)
        self.status_lbl = tk.Label(self.status_bar,
                                       text=self.T['welcome_msg'],
                                       bg=biosim_theme('panel'), fg=biosim_theme('text_dim'),
                                       font=('Arial', 9), anchor='w', padx=10)
        self.status_lbl.pack(side='left')
        self.seq_info_lbl = tk.Label(self.status_bar, text='',
                                          bg=biosim_theme('panel'), fg='#0550ae',
                                          font=('Courier', 9), padx=10)
        self.seq_info_lbl.pack(side='right')

    # ─── Tab builders (placeholders for Mol-0, real impl in Mol-1..5) ───

    def _build_tab_sequence(self):
        """Mol-1: Full Sequence Editor with ORF, translation, GC analysis."""
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=self.T['tab_seq'])
        self.tab_sequence = tab

        # ─── Top toolbar ───
        toolbar = tk.Frame(tab, bg=biosim_theme('panel'), height=42)
        toolbar.pack(fill='x', padx=4, pady=4)
        toolbar.pack_propagate(False)

        # Name input
        tk.Label(toolbar, text=('Имя:' if self.lang == 'ru' else 'Name:'),
                  bg=biosim_theme('panel'), fg=biosim_theme('text_dim'),
                  font=('Arial', 10)).pack(side='left', padx=(8, 2), pady=8)
        self.seq_name_var = tk.StringVar(value=self.sequence_name)
        tk.Entry(toolbar, textvariable=self.seq_name_var,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 10), width=20,
                  highlightthickness=1, highlightbackground=biosim_theme('border')
                  ).pack(side='left', padx=2, pady=8)

        # Type selector
        tk.Label(toolbar, text=('Тип:' if self.lang == 'ru' else 'Type:'),
                  bg=biosim_theme('panel'), fg=biosim_theme('text_dim'),
                  font=('Arial', 10)).pack(side='left', padx=(10, 2), pady=8)
        self.seq_type_var = tk.StringVar(value=self.sequence_type)
        type_dd = ttk.Combobox(toolbar, textvariable=self.seq_type_var,
                                  values=['dna', 'rna', 'protein'],
                                  state='readonly', width=8,
                                  font=('Arial', 10))
        type_dd.pack(side='left', padx=2, pady=8)
        type_dd.bind('<<ComboboxSelected>>', lambda e: self._on_seq_type_change())

        # Action buttons
        for icon, txt, color, cmd in [
            ('', 'RevComp', '#8250df', self._seq_revcomp),
            ('', 'Translate', '#1a7f37', self._seq_translate),
            ('', 'Find ORFs', '#0969da', self._seq_find_orfs),
            ('', 'Copy', '#9a6700', self._seq_copy),
            ('', 'Save', '#0969da', self._seq_save),
        ]:
            tk.Button(toolbar, text=icon + ' ' + txt,
                       bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=color, relief='solid', bd=1,
                       font=('Arial', 9, 'bold'), cursor='hand2',
                       padx=8, pady=2,
                       command=cmd).pack(side='left', padx=2, pady=6)

        # ─── Main split: editor + info ───
        main = tk.Frame(tab, bg=biosim_theme('bg'))
        main.pack(fill='both', expand=True, padx=4, pady=2)

        # Left: sequence editor
        left_frame = tk.LabelFrame(main,
                                       text='  ' + ('Последовательность' if self.lang == 'ru'
                                                        else 'Sequence') + '  ',
                                       bg=biosim_theme('bg'), fg='#0550ae',
                                       font=('Arial', 10, 'bold'), relief='flat')
        left_frame.pack(side='left', fill='both', expand=True, padx=(0, 4))

        editor_frame = tk.Frame(left_frame, bg=biosim_theme('panel'))
        editor_frame.pack(fill='both', expand=True, padx=4, pady=4)
        editor_sb = tk.Scrollbar(editor_frame)
        editor_sb.pack(side='right', fill='y')
        self.seq_text = tk.Text(editor_frame, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                   insertbackground='#1a7f37',
                                   font=('Courier', 11), relief='flat',
                                   wrap='none', highlightthickness=0,
                                   padx=8, pady=8,
                                   yscrollcommand=editor_sb.set)
        self.seq_text.pack(side='left', fill='both', expand=True)
        editor_sb.config(command=self.seq_text.yview)

        # Configure tags for highlighting
        self.seq_text.tag_config('A', foreground='#1a7f37')   # green
        self.seq_text.tag_config('T', foreground='#cf222e')   # red
        self.seq_text.tag_config('G', foreground='#9a6700')   # yellow
        self.seq_text.tag_config('C', foreground='#0969da')   # blue
        self.seq_text.tag_config('U', foreground='#cf222e')   # red (like T)
        self.seq_text.tag_config('N', foreground=biosim_theme('text_dim'))   # gray
        self.seq_text.tag_config('start_codon', background='#0969da', foreground='white')
        self.seq_text.tag_config('stop_codon', background='#cf222e', foreground='white')
        self.seq_text.tag_config('restriction', background='#1a7f37', foreground='#f5f7fa')
        self.seq_text.tag_config('orf_highlight', background='#eef1f5')

        # Bind text changes
        self.seq_text.bind('<KeyRelease>', self._on_seq_text_change)

        # Load current sequence if any
        if self.current_sequence:
            self.seq_text.insert('1.0', self.current_sequence)
            self._refresh_seq_highlights()

        # Right: info panel
        right_frame = tk.Frame(main, bg=biosim_theme('bg'), width=320)
        right_frame.pack(side='right', fill='y')
        right_frame.pack_propagate(False)

        # ── Cross-section transfer hub: MolBio ↔ SynBio sections ──
        link_frame = tk.Frame(right_frame, bg=biosim_theme('bg'))
        link_frame.pack(fill='x', pady=(0, 4))
        hub_btn = tk.Button(link_frame,
                  text=('🔗 Перенос / Синхронизация ▾' if self.lang == 'ru'
                        else '🔗 Transfer / Sync ▾'),
                  bg='#8250df', highlightbackground='#8250df', fg='white',
                  relief='flat', font=('Arial', 9, 'bold'), cursor='hand2',
                  pady=4)

        def show_molbio_hub():
            is_ru = (self.lang == 'ru')
            menu = tk.Menu(self, tearoff=0)
            menu.add_command(label=('— Перенести в: —' if is_ru else '— Send to: —'),
                             state='disabled')
            menu.add_command(label=('   CRISPR (найти мишени)' if is_ru
                                    else '   CRISPR (find targets)'),
                             command=self._molbio_link_to_crispr)
            menu.add_command(label=('   Конструктор (как часть)' if is_ru
                                    else '   Constructor (as part)'),
                             command=self._molbio_link_to_constructor)
            try:
                x = hub_btn.winfo_rootx()
                y = hub_btn.winfo_rooty() + hub_btn.winfo_height()
                menu.tk_popup(x, y)
            finally:
                menu.grab_release()

        hub_btn.config(command=show_molbio_hub)
        hub_btn.pack(fill='x')

        # Sanger chromatogram — opens a viewer for the current sequence
        tk.Button(right_frame,
                  text=('📊 Хроматограмма Сэнгера' if self.lang == 'ru'
                        else '📊 Sanger chromatogram'),
                  bg=biosim_theme('panel'), highlightbackground='#0969da',
                  fg='#0969da', relief='solid', bd=1,
                  font=('Arial', 9, 'bold'), cursor='hand2',
                  command=self._open_sanger_dialog).pack(fill='x', pady=(4, 0))

        # Sequence stats
        stats_frame = tk.LabelFrame(right_frame,
                                       text='  ' + ('Статистика' if self.lang == 'ru'
                                                        else 'Statistics') + '  ',
                                       bg=biosim_theme('bg'), fg='#1a7f37',
                                       font=('Arial', 10, 'bold'), relief='flat')
        stats_frame.pack(fill='x', pady=2)
        self.seq_stats_lbl = tk.Label(stats_frame, text='No sequence',
                                          bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                          font=('Courier', 9), anchor='w',
                                          justify='left')
        self.seq_stats_lbl.pack(padx=8, pady=4, anchor='w', fill='x')

        # Composition (A/T/G/C %)
        comp_frame = tk.LabelFrame(right_frame,
                                      text='  ' + ('Состав' if self.lang == 'ru'
                                                       else 'Composition') + '  ',
                                      bg=biosim_theme('bg'), fg='#0969da',
                                      font=('Arial', 10, 'bold'), relief='flat')
        comp_frame.pack(fill='x', pady=2)
        self.seq_comp_canvas = tk.Canvas(comp_frame, height=60,
                                            bg=biosim_theme('panel'), highlightthickness=0)
        self.seq_comp_canvas.pack(fill='x', padx=4, pady=4)
        self.seq_comp_lbl = tk.Label(comp_frame, text='',
                                         bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                         font=('Courier', 9), anchor='w',
                                         justify='left')
        self.seq_comp_lbl.pack(padx=8, pady=2, anchor='w')

        # ORFs list
        orfs_frame = tk.LabelFrame(right_frame,
                                       text='  ' + ('Открытые рамки' if self.lang == 'ru'
                                                        else 'ORFs') + '  ',
                                       bg=biosim_theme('bg'), fg='#8250df',
                                       font=('Arial', 10, 'bold'), relief='flat')
        orfs_frame.pack(fill='both', expand=True, pady=2)
        orfs_sb = tk.Scrollbar(orfs_frame)
        orfs_sb.pack(side='right', fill='y')
        self.orfs_listbox = tk.Listbox(orfs_frame,
                                          bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                          font=('Courier', 9), relief='flat',
                                          highlightthickness=0,
                                          selectbackground='#0969da',
                                          yscrollcommand=orfs_sb.set,
                                          height=12)
        self.orfs_listbox.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        orfs_sb.config(command=self.orfs_listbox.yview)
        self.orfs_listbox.bind('<<ListboxSelect>>', self._on_orf_select)
        self.orfs_data = []  # parallel to listbox

        # Bottom: translation panel (collapsible)
        trans_frame = tk.LabelFrame(tab,
                                       text='  ' + ('Трансляция (6 рамок)' if self.lang == 'ru'
                                                        else 'Translation (6-frame)') + '  ',
                                       bg=biosim_theme('bg'), fg='#9a6700',
                                       font=('Arial', 10, 'bold'), relief='flat')
        trans_frame.pack(fill='x', padx=4, pady=4)
        trans_inner = tk.Frame(trans_frame, bg=biosim_theme('panel'))
        trans_inner.pack(fill='x', padx=2, pady=2)
        self.trans_text = tk.Text(trans_inner, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                      font=('Courier', 10), relief='flat',
                                      wrap='none', highlightthickness=0,
                                      height=8, padx=8, pady=4)
        self.trans_text.pack(side='left', fill='x', expand=True)
        trans_sb = tk.Scrollbar(trans_inner, orient='horizontal')
        trans_sb.pack(side='bottom', fill='x')
        self.trans_text.config(xscrollcommand=trans_sb.set)
        trans_sb.config(command=self.trans_text.xview)
        self.trans_text.tag_config('frame_label', foreground='#0550ae',
                                       font=('Arial', 9, 'bold'))
        self.trans_text.tag_config('start_aa', background='#1a7f37', foreground='white')
        self.trans_text.tag_config('stop_aa', background='#cf222e', foreground='white')

        # Initial state
        self._refresh_seq_analysis()

    # ─── Sequence Editor methods ───

    def _molbio_link_to_crispr(self):
        """Send the current sequence to SynBio's CRISPR section and switch to it.
        Closes focus to SynBio window, opens it on the CRISPR tab with the
        sequence loaded — ready to scan for guides."""
        is_ru = (self.lang == 'ru')
        seq = self.seq_text.get('1.0', 'end').strip().upper()
        seq = ''.join(c for c in seq if c in 'ATGCUN')
        if len(seq) < 20:
            from tkinter import messagebox
            messagebox.showinfo('VIREON',
                ('Нужна последовательность минимум 20 нуклеотидов.' if is_ru
                 else 'Need a sequence of at least 20 nucleotides.'))
            return
        app = self.parent_main
        try:
            # Open (or focus) the SynBio window
            if hasattr(app, '_open_synbio'):
                app._open_synbio()
            sw = getattr(app, '_synbio_window', None)
            if sw is not None:
                # Switch to CRISPR tab and load the sequence
                sw._switch_tab('crispr')
                if hasattr(sw, 'load_external_crispr_sequence'):
                    sw.load_external_crispr_sequence(seq)
                sw.lift(); sw.focus_force()
                # Optionally lower this MolBio window
                self.lower()
        except Exception as e:
            from tkinter import messagebox
            messagebox.showerror('VIREON', 'SynBio: ' + str(e)[:80])

    def load_external_sequence(self, sequence, target_tab=None):
        """Receive a sequence from another section (e.g. SynBio CRISPR) and
        load it into the MolBio sequence editor."""
        try:
            seq = ''.join(c for c in (sequence or '').upper() if c in 'ATGCUN')
            if hasattr(self, 'seq_text') and seq:
                self.seq_text.delete('1.0', 'end')
                self.seq_text.insert('1.0', seq)
                self.current_sequence = seq
                self._refresh_seq_highlights()
                self._refresh_seq_analysis()
                # Switch to the sequence tab if possible
                try:
                    self.notebook.select(0)
                except Exception as _exc:
                    biosim_log('warn', 'load_external_sequence: {}'.format(_exc))
        except Exception as _exc:
            biosim_log('warn', 'load_external_sequence: {}'.format(_exc))

    def get_current_sequence(self):
        """Return the current (possibly edited) sequence — used by other
        sections to sync changes back (bidirectional transfer)."""
        try:
            if hasattr(self, 'seq_text'):
                raw = self.seq_text.get('1.0', 'end').strip()
                return ''.join(c for c in raw.upper() if c in 'ATGCUN')
            return getattr(self, 'current_sequence', '') or ''
        except Exception as _exc:
            biosim_log('warn', 'get_current_sequence: {}'.format(_exc))
            return ''

    def _molbio_link_to_constructor(self):
        """Send current sequence to SynBio Constructor (as a gene part)."""
        is_ru = (self.lang == 'ru')
        seq = self.get_current_sequence()
        if len(seq) < 10:
            from tkinter import messagebox
            messagebox.showinfo('VIREON',
                ('Нужна последовательность ≥10 нт.' if is_ru
                 else 'Need a sequence of ≥10 nt.'))
            return
        app = self.parent_main
        try:
            if hasattr(app, '_open_synbio'):
                app._open_synbio()
            sw = getattr(app, '_synbio_window', None)
            if sw is not None:
                sw._switch_tab('constructor')
                if hasattr(sw, 'load_external_constructor_sequence'):
                    sw.load_external_constructor_sequence(seq)
                sw.lift(); sw.focus_force()
                self.lower()
        except Exception as e:
            from tkinter import messagebox
            messagebox.showerror('VIREON', 'SynBio: ' + str(e)[:80])

    def _on_seq_text_change(self, event=None):
        """Called when user types — auto-clean and update analysis."""
        try:
            raw = self.seq_text.get('1.0', 'end-1c')
            # Clean: keep only letters
            cleaned = ''.join(c for c in raw if c.isalpha()).upper()
            if cleaned != raw.replace('\n', '').replace(' ', '').upper():
                # Re-insert cleaned
                cursor_pos = self.seq_text.index('insert')
                self.seq_text.delete('1.0', 'end')
                self.seq_text.insert('1.0', cleaned)
                try:
                    self.seq_text.mark_set('insert', cursor_pos)
                except Exception as _exc:
                    biosim_log('warn', '_on_seq_text_change: {}'.format(_exc))
            self.current_sequence = cleaned
            # Auto-detect type
            detected = detect_sequence_type(cleaned)
            if detected != self.sequence_type:
                self.sequence_type = detected
                self.seq_type_var.set(detected)
            self._refresh_seq_analysis()
            self._refresh_seq_highlights()
        except Exception as e:
            biosim_log('error', 'seq_text_change: ' + str(e))

    def _on_seq_type_change(self):
        """User manually changed type."""
        self.sequence_type = self.seq_type_var.get()
        self._refresh_seq_analysis()

    def _refresh_seq_highlights(self):
        """Add colored highlighting to nucleotides."""
        try:
            # Remove old tags
            for tag in ('A', 'T', 'G', 'C', 'U', 'N',
                          'start_codon', 'stop_codon'):
                self.seq_text.tag_remove(tag, '1.0', 'end')
            seq = self.current_sequence
            if not seq:
                return
            # Tag each character (limit to first 2000 for performance)
            limit = min(len(seq), 2000)
            for i, c in enumerate(seq[:limit]):
                if c in 'ATGCUN':
                    start = '1.0+{}c'.format(i)
                    end = '1.0+{}c'.format(i + 1)
                    self.seq_text.tag_add(c, start, end)
        except Exception as e:
            biosim_log('error', 'refresh_highlights: ' + str(e))

    def _refresh_seq_analysis(self):
        """Update stats, composition, ORFs."""
        seq = self.current_sequence
        if not seq:
            self.seq_stats_lbl.config(text='No sequence')
            self.seq_comp_lbl.config(text='')
            self.seq_comp_canvas.delete('all')
            self.orfs_listbox.delete(0, 'end')
            self.trans_text.delete('1.0', 'end')
            return

        # Stats
        stats = []
        stats.append('Length: {} bp'.format(len(seq)) if self.sequence_type != 'protein'
                       else 'Length: {} aa'.format(len(seq)))

        if self.sequence_type in ('dna', 'rna'):
            gc = gc_content(seq)
            stats.append('GC: {:.2f}%'.format(gc))
            stats.append('AT: {:.2f}%'.format(100 - gc))
            stats.append('MW (dsDNA): {:.0f} Da'.format(molecular_weight_dna(seq)))
            stats.append('Tm (estimate): {:.1f}°C'.format(
                self._estimate_tm(seq)))
        elif self.sequence_type == 'protein':
            mw = molecular_weight_protein(seq)
            stats.append('MW: {:.0f} Da'.format(mw))
            stats.append('Avg MW: {:.1f}'.format(mw / len(seq) if seq else 0))
            # Charge at pH 7
            charge = sum(AA_PROPERTIES.get(aa, {'charge': 0})['charge']
                          for aa in seq)
            stats.append('Net charge (pH 7): {:.1f}'.format(charge))

        self.seq_stats_lbl.config(text='\n'.join(stats))

        # Composition bar chart
        self._draw_composition_chart()

        # Translation (only if DNA/RNA)
        if self.sequence_type in ('dna', 'rna'):
            self._show_translation_6frame()
        else:
            self.trans_text.delete('1.0', 'end')

        # Update top-level status
        self._update_status()

    def _estimate_tm(self, seq):
        """Simple Wallace rule Tm estimate."""
        if not seq:
            return 0.0
        seq = seq.upper()
        if len(seq) < 14:
            # Wallace: 2*(A+T) + 4*(G+C)
            at = seq.count('A') + seq.count('T') + seq.count('U')
            gc = seq.count('G') + seq.count('C')
            return 2 * at + 4 * gc
        # Long: GC method
        gc = gc_content(seq)
        return 64.9 + 41 * (gc / 100 - 0.1644)

    def _draw_composition_chart(self):
        """Draw nucleotide composition as bars."""
        try:
            self.seq_comp_canvas.delete('all')
            seq = self.current_sequence
            if not seq:
                return
            seq_u = seq.upper()
            total = len(seq_u)
            if total == 0:
                return

            counts = {}
            colors = {'A': '#1a7f37', 'T': '#cf222e', 'G': '#9a6700',
                       'C': '#0969da', 'U': '#cf222e'}
            for ch in 'ATGCU':
                counts[ch] = seq_u.count(ch)

            # Draw bars
            w = int(self.seq_comp_canvas.winfo_width()) or 280
            chars = [ch for ch in 'ATGCU' if counts[ch] > 0]
            if not chars:
                return
            bar_w = (w - 20) / len(chars)
            comp_lines = []
            for i, ch in enumerate(chars):
                pct = 100 * counts[ch] / total
                bar_h = int(40 * pct / 100)
                x = 10 + i * bar_w
                self.seq_comp_canvas.create_rectangle(
                    x, 50 - bar_h, x + bar_w - 4, 50,
                    fill=colors.get(ch, '#57606a'), outline='')
                self.seq_comp_canvas.create_text(
                    x + bar_w / 2, 55,
                    text=ch, fill='#1c2128',
                    font=('Arial', 9, 'bold'))
                self.seq_comp_canvas.create_text(
                    x + bar_w / 2, 50 - bar_h - 6,
                    text='{:.0f}%'.format(pct),
                    fill='#57606a', font=('Arial', 8))
                comp_lines.append('{}: {} ({:.1f}%)'.format(ch, counts[ch], pct))
            self.seq_comp_lbl.config(text='   '.join(comp_lines))
        except Exception as e:
            biosim_log('error', 'comp chart: ' + str(e))

    def _show_translation_6frame(self):
        """Display 6-frame translation in trans_text."""
        try:
            self.trans_text.delete('1.0', 'end')
            seq = self.current_sequence
            if not seq or len(seq) < 3:
                return
            frames = translate_six_frames(seq)
            labels = ['+1', '+2', '+3', '-1', '-2', '-3']
            for i, (lbl, prot) in enumerate(zip(labels, frames)):
                self.trans_text.insert('end', '{}: '.format(lbl), 'frame_label')
                # Insert with highlighting
                for ch in prot[:200]:  # limit display
                    if ch == 'M':
                        self.trans_text.insert('end', ch, 'start_aa')
                    elif ch == '*':
                        self.trans_text.insert('end', ch, 'stop_aa')
                    else:
                        self.trans_text.insert('end', ch)
                if len(prot) > 200:
                    self.trans_text.insert('end', '... (' + str(len(prot) - 200) + ' more)')
                self.trans_text.insert('end', '\n')
        except Exception as e:
            biosim_log('error', 'show_translation: ' + str(e))

    def _seq_revcomp(self):
        """Reverse complement current sequence."""
        if not self.current_sequence:
            return
        if self.sequence_type != 'dna':
            self.status_lbl.config(text='Reverse complement only for DNA',
                                       fg='#9a6700')
            return
        rc = reverse_complement(self.current_sequence)
        self.seq_text.delete('1.0', 'end')
        self.seq_text.insert('1.0', rc)
        self.current_sequence = rc
        self._refresh_seq_analysis()
        self._refresh_seq_highlights()
        self.status_lbl.config(text='Reverse complemented', fg='#1a7f37')

    def _seq_translate(self):
        """Replace current sequence with its protein translation."""
        if not self.current_sequence:
            return
        if self.sequence_type not in ('dna', 'rna'):
            self.status_lbl.config(text='Translation requires DNA/RNA',
                                       fg='#9a6700')
            return
        protein = translate_dna(self.current_sequence, 0).rstrip('*')
        # Open new dialog showing the protein
        dlg = tk.Toplevel(self)
        dlg.title('Translation result')
        dlg.geometry('700x500')
        dlg.configure(bg=biosim_theme('bg'))
        dlg.transient(self)
        tk.Label(dlg, text='Translation (frame +1)',
                  bg=biosim_theme('bg'), fg='#1a7f37',
                  font=('Arial', 14, 'bold')).pack(pady=8)
        tk.Label(dlg, text='{} aa · MW {:.0f} Da'.format(
                    len(protein), molecular_weight_protein(protein)),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 10, 'italic')).pack()
        prot_txt = tk.Text(dlg, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                              font=('Courier', 11), relief='flat',
                              wrap='word', padx=10, pady=10)
        prot_txt.pack(fill='both', expand=True, padx=14, pady=8)
        # Insert with color by AA type
        for ch in protein:
            tag = AA_PROPERTIES.get(ch, {}).get('type', '')
            color = AA_COLORS.get(tag, '#1c2128')
            tag_name = 'aa_' + tag
            try:
                prot_txt.tag_config(tag_name, foreground=color)
            except Exception as _exc:
                biosim_log('warn', '_seq_translate: {}'.format(_exc))
            prot_txt.insert('end', ch, tag_name)
        prot_txt.config(state='disabled')

        def copy_protein():
            self.clipboard_clear()
            self.clipboard_append(protein)
            self.status_lbl.config(text='Protein copied', fg='#1a7f37')

        def load_as_protein():
            self.current_sequence = protein
            self.sequence_type = 'protein'
            self.seq_type_var.set('protein')
            self.seq_text.delete('1.0', 'end')
            self.seq_text.insert('1.0', protein)
            self._refresh_seq_analysis()
            dlg.destroy()

        btn_row = tk.Frame(dlg, bg=biosim_theme('bg'))
        btn_row.pack(fill='x', padx=14, pady=8)
        tk.Button(btn_row, text='Copy',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=12, pady=4, command=copy_protein).pack(side='left', padx=2)
        tk.Button(btn_row, text='Load as protein',
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=12, pady=4, command=load_as_protein).pack(side='left', padx=2)
        tk.Button(btn_row, text='Close',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2',
                   padx=12, pady=4, command=dlg.destroy).pack(side='right', padx=2)

    def _seq_find_orfs(self):
        """Find and list all ORFs."""
        if not self.current_sequence:
            self.status_lbl.config(text='No sequence to analyze', fg='#9a6700')
            return
        if self.sequence_type not in ('dna', 'rna'):
            self.status_lbl.config(text='ORF finder requires DNA/RNA',
                                       fg='#9a6700')
            return
        # Ask for minimum length
        from tkinter import simpledialog
        try:
            min_aa = simpledialog.askinteger(
                'Find ORFs', 'Minimum ORF length (amino acids):',
                initialvalue=30, minvalue=1, maxvalue=10000,
                parent=self)
        except Exception as _exc:
            biosim_log('warn', '_seq_find_orfs: {}'.format(_exc))
            min_aa = 30
        if min_aa is None:
            return

        orfs = find_orfs(self.current_sequence, min_length=min_aa)
        self.orfs_data = orfs
        self.orfs_listbox.delete(0, 'end')
        if not orfs:
            self.orfs_listbox.insert('end', '(no ORFs found)')
            self.status_lbl.config(text='No ORFs ≥ {} aa'.format(min_aa),
                                       fg='#9a6700')
            return
        for i, orf in enumerate(orfs):
            txt = '{}{}: {}-{} ({} aa)'.format(
                orf['strand'], abs(orf['frame']),
                orf['start'] + 1, orf['end'],
                orf['length_aa'])
            self.orfs_listbox.insert('end', txt)
        self.status_lbl.config(text='Found ' + str(len(orfs)) + ' ORFs',
                                   fg='#1a7f37')

    def _on_orf_select(self, event=None):
        """User clicked on an ORF — show details."""
        sel = self.orfs_listbox.curselection()
        if not sel or sel[0] >= len(self.orfs_data):
            return
        orf = self.orfs_data[sel[0]]
        # Show details dialog
        dlg = tk.Toplevel(self)
        dlg.title('ORF Details')
        dlg.geometry('640x500')
        dlg.configure(bg=biosim_theme('bg'))
        dlg.transient(self)

        tk.Label(dlg, text='ORF Details',
                  bg=biosim_theme('bg'), fg='#8250df',
                  font=('Arial', 14, 'bold')).pack(pady=8)

        info = 'Strand: {}\n'.format(orf['strand'])
        info += 'Frame: {}\n'.format(orf['frame'])
        info += 'Start: {} (nt)\n'.format(orf['start'] + 1)
        info += 'End: {} (nt)\n'.format(orf['end'])
        info += 'Length: {} nt = {} aa\n'.format(orf['length_nt'], orf['length_aa'])
        info += 'Protein MW: {:.0f} Da'.format(molecular_weight_protein(orf['protein']))
        tk.Label(dlg, text=info, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                  font=('Courier', 10), justify='left').pack(padx=14, pady=4, anchor='w')

        tk.Label(dlg, text='Protein:', bg=biosim_theme('bg'), fg='#1a7f37',
                  font=('Arial', 11, 'bold')).pack(anchor='w', padx=14, pady=(8, 2))
        prot_txt = tk.Text(dlg, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                              font=('Courier', 10), relief='flat',
                              wrap='word', height=8, padx=10, pady=8)
        prot_txt.pack(fill='x', padx=14, pady=4)
        prot_txt.insert('1.0', orf['protein'])
        prot_txt.config(state='disabled')

        tk.Button(dlg, text='Copy protein',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=lambda: (self.clipboard_clear(),
                                       self.clipboard_append(orf['protein']))
                   ).pack(fill='x', padx=14, pady=4)
        tk.Button(dlg, text='Close',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy).pack(fill='x', padx=14, pady=(4, 14))

    def _seq_copy(self):
        """Copy current sequence to clipboard."""
        if not self.current_sequence:
            return
        self.clipboard_clear()
        self.clipboard_append(self.current_sequence)
        self.status_lbl.config(text='Sequence copied', fg='#1a7f37')

    def _seq_save(self):
        """Save sequence as FASTA file."""
        if not self.current_sequence:
            return
        from tkinter import filedialog
        filepath = filedialog.asksaveasfilename(
            title='Save FASTA',
            defaultextension='.fasta',
            initialfile=self.seq_name_var.get() + '.fasta',
            filetypes=[('FASTA', '*.fasta *.fa'), ('All files', '*.*')])
        if not filepath:
            return
        try:
            name = self.seq_name_var.get() or 'Untitled'
            content = '>' + name + '\n' + format_sequence_blocks(
                self.current_sequence, block_size=10, line_size=60,
                show_position=False) + '\n'
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
            self.status_lbl.config(text='Saved · ' + os.path.basename(filepath),
                                       fg='#1a7f37')
        except Exception as e:
            biosim_log('warn', '_seq_save: {}'.format(e))
            self.status_lbl.config(text='Save error: ' + str(e)[:50],
                                       fg='#cf222e')

    def _build_tab_primer(self):
        """Mol-2: Full Primer Designer."""
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=self.T['tab_primer'])
        self.tab_primer = tab

        # ─── Top: controls ───
        ctrl = tk.LabelFrame(tab, text='  ' + ('Дизайн праймеров' if self.lang == 'ru'
                                                       else 'Primer Design') + '  ',
                                bg=biosim_theme('bg'), fg='#8250df',
                                font=('Arial', 11, 'bold'), relief='flat')
        ctrl.pack(fill='x', padx=4, pady=4)

        # Mode buttons
        mode_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        mode_frame.pack(fill='x', padx=8, pady=6)
        tk.Label(mode_frame, text=('Режим:' if self.lang == 'ru' else 'Mode:'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 10, 'bold')).pack(side='left', padx=(0, 8))

        self.primer_mode_var = tk.StringVar(value='design')
        for val, lbl, col in [
            ('design', '' + ('Авто-дизайн' if self.lang == 'ru' else 'Auto-design'), '#1a7f37'),
            ('analyze', '' + ('Анализ праймера' if self.lang == 'ru' else 'Analyze primer'), '#0969da'),
            ('pair', '' + ('Анализ пары' if self.lang == 'ru' else 'Analyze pair'), '#9a6700'),
        ]:
            tk.Radiobutton(mode_frame, text=lbl, value=val,
                            variable=self.primer_mode_var,
                            bg=biosim_theme('bg'), fg=col,
                            selectcolor='#eef1f5',
                            activebackground='#f5f7fa',
                            font=('Arial', 9, 'bold'),
                            command=self._switch_primer_mode
                            ).pack(side='left', padx=4)

        # Parameter grid
        params_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        params_frame.pack(fill='x', padx=8, pady=4)

        # Row 1: target range (for design mode)
        tk.Label(params_frame, text=('Целевая область (start-end):'
                                        if self.lang == 'ru' else
                                        'Target region (start-end):'),
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 9)).grid(row=0, column=0, sticky='w', padx=4)
        self.primer_target_start = tk.IntVar(value=1)
        self.primer_target_end = tk.IntVar(value=500)
        tk.Entry(params_frame, textvariable=self.primer_target_start,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=8
                  ).grid(row=0, column=1, padx=2)
        tk.Label(params_frame, text=' — ', bg=biosim_theme('bg'), fg=biosim_theme('text_dim')
                  ).grid(row=0, column=2)
        tk.Entry(params_frame, textvariable=self.primer_target_end,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=8
                  ).grid(row=0, column=3, padx=2)

        # Row 2: Tm range
        tk.Label(params_frame, text='Tm (°C):',
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 9)).grid(row=1, column=0, sticky='w', padx=4)
        self.primer_tm_min = tk.IntVar(value=55)
        self.primer_tm_max = tk.IntVar(value=65)
        tk.Entry(params_frame, textvariable=self.primer_tm_min,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=8
                  ).grid(row=1, column=1, padx=2)
        tk.Label(params_frame, text=' — ', bg=biosim_theme('bg'), fg=biosim_theme('text_dim')
                  ).grid(row=1, column=2)
        tk.Entry(params_frame, textvariable=self.primer_tm_max,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=8
                  ).grid(row=1, column=3, padx=2)

        # Row 3: GC range
        tk.Label(params_frame, text='GC (%):',
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 9)).grid(row=2, column=0, sticky='w', padx=4)
        self.primer_gc_min = tk.IntVar(value=40)
        self.primer_gc_max = tk.IntVar(value=60)
        tk.Entry(params_frame, textvariable=self.primer_gc_min,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=8
                  ).grid(row=2, column=1, padx=2)
        tk.Label(params_frame, text=' — ', bg=biosim_theme('bg'), fg=biosim_theme('text_dim')
                  ).grid(row=2, column=2)
        tk.Entry(params_frame, textvariable=self.primer_gc_max,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=8
                  ).grid(row=2, column=3, padx=2)

        # Row 4: length
        tk.Label(params_frame, text=('Длина (nt):' if self.lang == 'ru'
                                       else 'Length (nt):'),
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 9)).grid(row=3, column=0, sticky='w', padx=4)
        self.primer_len_min = tk.IntVar(value=18)
        self.primer_len_max = tk.IntVar(value=25)
        tk.Entry(params_frame, textvariable=self.primer_len_min,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=8
                  ).grid(row=3, column=1, padx=2)
        tk.Label(params_frame, text=' — ', bg=biosim_theme('bg'), fg=biosim_theme('text_dim')
                  ).grid(row=3, column=2)
        tk.Entry(params_frame, textvariable=self.primer_len_max,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=8
                  ).grid(row=3, column=3, padx=2)

        # Run button
        tk.Button(params_frame, text='▶ ' + ('Найти праймеры' if self.lang == 'ru'
                                                 else 'Find primers'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4,
                   command=self._run_primer_design
                   ).grid(row=0, column=5, rowspan=2, padx=12, sticky='e')

        # Manual primer input (for analyze mode)
        self.primer_manual_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        # Hidden initially — shown in analyze mode
        tk.Label(self.primer_manual_frame, text='Forward (5\'3\'):',
                  bg=biosim_theme('bg'), fg='#1a7f37',
                  font=('Arial', 9, 'bold')).pack(anchor='w', padx=4)
        self.primer_fw_var = tk.StringVar()
        tk.Entry(self.primer_manual_frame, textvariable=self.primer_fw_var,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Courier', 10),
                  highlightthickness=1, highlightbackground=biosim_theme('border')
                  ).pack(fill='x', padx=4, pady=(0, 4), ipady=3)
        tk.Label(self.primer_manual_frame, text='Reverse (5\'3\'):',
                  bg=biosim_theme('bg'), fg='#cf222e',
                  font=('Arial', 9, 'bold')).pack(anchor='w', padx=4)
        self.primer_rv_var = tk.StringVar()
        tk.Entry(self.primer_manual_frame, textvariable=self.primer_rv_var,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Courier', 10),
                  highlightthickness=1, highlightbackground=biosim_theme('border')
                  ).pack(fill='x', padx=4, pady=(0, 4), ipady=3)

        # ─── Results area ───
        results_frame = tk.LabelFrame(tab,
                                          text='  ' + ('Результаты' if self.lang == 'ru'
                                                           else 'Results') + '  ',
                                          bg=biosim_theme('bg'), fg='#1a7f37',
                                          font=('Arial', 11, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=4, pady=4)

        # Scrollable results
        canv = tk.Canvas(results_frame, bg=biosim_theme('bg'), highlightthickness=0)
        sb = tk.Scrollbar(results_frame, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        sb.pack(side='right', fill='y')
        self.primer_results_inner = tk.Frame(canv, bg=biosim_theme('bg'))
        self.primer_results_inner.bind('<Configure>',
                                            lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=self.primer_results_inner,
                             anchor='nw', width=1100)

        # Initial message
        self._primer_show_initial_msg()

    def _switch_primer_mode(self):
        """Toggle between design / analyze / pair modes."""
        mode = self.primer_mode_var.get()
        if mode in ('analyze', 'pair'):
            self.primer_manual_frame.pack(fill='x', padx=8, pady=4, after=self.primer_manual_frame.master.winfo_children()[2])
        else:
            self.primer_manual_frame.pack_forget()
        self._primer_clear_results()
        self._primer_show_initial_msg()

    def _primer_clear_results(self):
        for w in self.primer_results_inner.winfo_children():
            w.destroy()

    def _primer_show_initial_msg(self):
        self._primer_clear_results()
        mode = self.primer_mode_var.get()
        msgs = {
            'design': ('Нажми "▶ Найти праймеры" для автодизайна'
                        if self.lang == 'ru' else
                        'Click "▶ Find primers" for auto-design'),
            'analyze': ('Введи forward праймер и жми ▶'
                          if self.lang == 'ru' else
                          'Enter forward primer and click ▶'),
            'pair': ('Введи forward и reverse праймеры и жми ▶'
                      if self.lang == 'ru' else
                      'Enter forward and reverse primers and click ▶'),
        }
        tk.Label(self.primer_results_inner, text=msgs.get(mode, ''),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 11, 'italic'),
                  padx=20, pady=20).pack(anchor='center')

    def _run_primer_design(self):
        """Run primer design / analysis based on current mode."""
        mode = self.primer_mode_var.get()
        if mode == 'design':
            self._do_primer_design()
        elif mode == 'analyze':
            self._do_primer_analyze()
        elif mode == 'pair':
            self._do_primer_pair()

    def _do_primer_design(self):
        """Auto-design primer pairs for target region."""
        if not self.current_sequence:
            self._primer_show_error('Load a DNA sequence first'
                                       if self.lang == 'en' else
                                       'Сначала загрузи ДНК последовательность')
            return
        if self.sequence_type not in ('dna', 'rna'):
            self._primer_show_error('Need DNA/RNA sequence'
                                       if self.lang == 'en' else
                                       'Нужна ДНК/РНК')
            return

        start = self.primer_target_start.get() - 1  # 1-indexed 0
        end = self.primer_target_end.get()
        if end <= start:
            self._primer_show_error('Invalid target range')
            return

        self._primer_clear_results()
        tk.Label(self.primer_results_inner, text='⏳ ' + ('Поиск праймеров...'
                                                              if self.lang == 'ru' else
                                                              'Searching for primers...'),
                  bg=biosim_theme('bg'), fg='#9a6700',
                  font=('Arial', 11, 'italic'), pady=10).pack()
        self.update()

        try:
            pairs = design_primer_pair(
                self.current_sequence, start, end,
                primer_length=(self.primer_len_min.get(), self.primer_len_max.get()),
                tm_range=(self.primer_tm_min.get(), self.primer_tm_max.get()),
                gc_range=(self.primer_gc_min.get(), self.primer_gc_max.get()),
            )
        except Exception as e:
            self._primer_show_error('Error: ' + str(e)[:80])
            biosim_log('error', 'primer design: ' + str(e))
            return

        self._primer_clear_results()
        if not pairs:
            self._primer_show_error(
                'No good primer pairs found. Try widening Tm or GC range.'
                if self.lang == 'en' else
                'Нет подходящих пар. Расширь диапазон Tm или GC.')
            return

        # Header
        tk.Label(self.primer_results_inner,
                  text=('Найдено пар' if self.lang == 'ru' else 'Found pairs') +
                        ': {}'.format(len(pairs)),
                  bg=biosim_theme('bg'), fg='#1a7f37',
                  font=('Arial', 12, 'bold'), pady=6).pack(anchor='w', padx=12)

        # Show top pairs as cards
        for i, pair in enumerate(pairs[:5], 1):
            self._render_primer_pair_card(pair, i)

    def _render_primer_pair_card(self, pair, idx):
        """Display one primer pair as a visual card."""
        card = tk.Frame(self.primer_results_inner, bg=biosim_theme('panel'),
                          highlightthickness=2,
                          highlightbackground=self._score_color(pair['combined_score']))
        card.pack(fill='x', padx=8, pady=4)

        # Header
        head = tk.Frame(card, bg=biosim_theme('panel'))
        head.pack(fill='x', padx=10, pady=(8, 4))
        tk.Label(head, text='#{} Pair'.format(idx),
                  bg=biosim_theme('panel'), fg='#0550ae',
                  font=('Arial', 12, 'bold')).pack(side='left')
        tk.Label(head, text=' · Score: {:.0f}'.format(pair['combined_score']),
                  bg=biosim_theme('panel'),
                  fg=self._score_color(pair['combined_score']),
                  font=('Arial', 11, 'bold')).pack(side='left', padx=4)
        tk.Label(head, text='Product: {} bp'.format(pair['product_length']),
                  bg=biosim_theme('panel'), fg='#9a6700',
                  font=('Courier', 10)).pack(side='right')
        tk.Label(head, text='ΔTm: {:.1f}°C'.format(pair['tm_diff']),
                  bg=biosim_theme('panel'), fg=biosim_theme('text_dim'),
                  font=('Courier', 9)).pack(side='right', padx=8)

        # Forward primer
        fw = pair['forward']
        self._render_single_primer_row(card, 'Forward', fw, '#1a7f37')

        # Reverse primer
        rv = pair['reverse']
        self._render_single_primer_row(card, 'Reverse', rv, '#cf222e')

        # Cross-dimer warning
        if pair.get('cross_dimers'):
            warn_text = '' + ('Cross-dimer: ' if self.lang == 'en' else 'Перекрёстный димер: ')
            for d in pair['cross_dimers'][:2]:
                warn_text += '{} bp ({}). '.format(d['match_length'],
                                                       '3-end' if d['is_3_end'] else 'mid')
            tk.Label(card, text=warn_text,
                      bg=biosim_theme('panel'), fg='#9a6700',
                      font=('Arial', 9, 'italic'),
                      wraplength=900, justify='left').pack(anchor='w', padx=12, pady=2)

        # Buttons
        btns = tk.Frame(card, bg=biosim_theme('panel'))
        btns.pack(fill='x', padx=10, pady=(2, 8))

        def copy_pair():
            self.clipboard_clear()
            self.clipboard_append('Forward: {}\nReverse: {}\nProduct: {} bp'.format(
                fw['sequence'], rv['sequence'], pair['product_length']))
            self.status_lbl.config(text='Primer pair copied', fg='#1a7f37')

        tk.Button(btns, text='Copy', bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'),
                   relief='solid', bd=1, font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=8, pady=2, command=copy_pair).pack(side='left', padx=2)

    def _render_single_primer_row(self, parent, label, primer_data, color):
        """Display one primer's properties in a row."""
        row = tk.Frame(parent, bg=biosim_theme('bg'))
        row.pack(fill='x', padx=10, pady=2)

        # Label
        tk.Label(row, text=label, bg=biosim_theme('bg'), fg=color,
                  font=('Arial', 10, 'bold'), width=12, anchor='w'
                  ).pack(side='left')

        # Sequence (monospace)
        seq_lbl = tk.Label(row, text=primer_data['sequence'],
                              bg=biosim_theme('bg'), fg=biosim_theme('text'),
                              font=('Courier', 11, 'bold'))
        seq_lbl.pack(side='left', padx=4)

        # Click to copy
        def copy_this(e):
            self.clipboard_clear()
            self.clipboard_append(primer_data['sequence'])
            self.status_lbl.config(text='' + label + ' copied', fg='#1a7f37')
        seq_lbl.bind('<Button-1>', copy_this)
        seq_lbl.config(cursor='hand2')

        # Properties
        props_text = ' · {} nt · Tm {:.1f}°C · GC {:.0f}%'.format(
            primer_data['length'], primer_data['tm'],
            primer_data['gc_content'])
        tk.Label(row, text=props_text,
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Courier', 9)).pack(side='left', padx=4)

        # Warnings
        if primer_data['hairpins']:
            tk.Label(row, text='HP', bg=biosim_theme('bg'), fg='#9a6700',
                      font=('Arial', 8, 'bold')).pack(side='left', padx=4)
        if primer_data.get('self_dimers'):
            sd = primer_data['self_dimers']
            if any(d['is_3_end'] for d in sd):
                tk.Label(row, text='SD', bg=biosim_theme('bg'), fg='#cf222e',
                          font=('Arial', 8, 'bold')).pack(side='left', padx=4)

    def _do_primer_analyze(self):
        """Analyze a single primer entered by user."""
        seq = ''.join(c.upper() for c in self.primer_fw_var.get() if c.isalpha())
        if not seq or len(seq) < 10:
            self._primer_show_error('Enter a primer (min 10 nt)')
            return
        if any(c not in 'ATGCU' for c in seq):
            self._primer_show_error('Invalid DNA characters')
            return

        result = score_primer(seq, tm_target=60)
        if not result:
            self._primer_show_error('Primer length must be 15-35 nt')
            return

        self._primer_clear_results()

        # Big info card
        card = tk.Frame(self.primer_results_inner, bg=biosim_theme('panel'),
                          highlightthickness=2,
                          highlightbackground=self._score_color(result['score']))
        card.pack(fill='x', padx=10, pady=10)

        # Header
        head = tk.Frame(card, bg=biosim_theme('panel'))
        head.pack(fill='x', padx=12, pady=(10, 4))
        tk.Label(head, text=('Анализ' if self.lang == 'ru' else 'Analysis'),
                  bg=biosim_theme('panel'), fg='#0550ae',
                  font=('Arial', 14, 'bold')).pack(side='left')
        tk.Label(head, text=' · Score: {:.0f}/100'.format(result['score']),
                  bg=biosim_theme('panel'),
                  fg=self._score_color(result['score']),
                  font=('Arial', 13, 'bold')).pack(side='left', padx=4)

        # Sequence
        tk.Label(card, text=seq, bg=biosim_theme('panel'), fg='#1a7f37',
                  font=('Courier', 13, 'bold')).pack(padx=12, pady=4)

        # Properties grid
        grid = tk.Frame(card, bg=biosim_theme('panel'))
        grid.pack(fill='x', padx=12, pady=4)

        rows = [
            ('Length:', '{} nt'.format(result['length']), '#1c2128'),
            ('Tm (NN):', '{:.1f} °C'.format(result['tm']),
              self._tm_color(result['tm'])),
            ('Tm (Wallace):', '{:.1f} °C'.format(result['tm_wallace']), '#57606a'),
            ('Tm (GC):', '{:.1f} °C'.format(result['tm_gc']), '#57606a'),
            ('GC content:', '{:.1f} %'.format(result['gc_content']),
              self._gc_color(result['gc_content'])),
            ('3\' clamp:', '{} GC in last {} ({})'.format(
                result['clamp']['gc_count'],
                len(result['clamp']['last_bases']),
                result['clamp']['status']),
              self._clamp_color(result['clamp']['status'])),
        ]
        for i, (k, v, c) in enumerate(rows):
            tk.Label(grid, text=k, bg=biosim_theme('panel'), fg=biosim_theme('text_dim'),
                      font=('Arial', 10), anchor='e', width=14
                      ).grid(row=i, column=0, sticky='e', padx=4, pady=1)
            tk.Label(grid, text=v, bg=biosim_theme('panel'), fg=c,
                      font=('Courier', 10, 'bold'), anchor='w'
                      ).grid(row=i, column=1, sticky='w', padx=4, pady=1)

        # Warnings section
        warn_frame = tk.LabelFrame(card,
                                       text='  ' + ('Предупреждения' if self.lang == 'ru'
                                                        else 'Warnings') + '  ',
                                       bg=biosim_theme('panel'), fg='#9a6700',
                                       font=('Arial', 9, 'bold'), relief='flat')
        warn_frame.pack(fill='x', padx=12, pady=4)

        any_warning = False
        if result['hairpins']:
            any_warning = True
            for h in result['hairpins'][:3]:
                tk.Label(warn_frame,
                          text='Hairpin: stem {} bp, loop {} bp, ΔG ≈ {:.1f} kcal/mol'.format(
                              h['stem_length'], h['loop_length'], h['dG_estimate']),
                          bg=biosim_theme('panel'), fg='#bc4c00',
                          font=('Courier', 9), anchor='w'
                          ).pack(anchor='w', padx=8, pady=1)
        if result['self_dimers']:
            for d in result['self_dimers'][:3]:
                any_warning = True
                marker = '3\'-END' if d['is_3_end'] else '~ internal'
                tk.Label(warn_frame,
                          text='Self-dimer: {} bp match ({}) — {}'.format(
                              d['match_length'], marker, d['matched_seq']),
                          bg=biosim_theme('panel'), fg='#cf222e',
                          font=('Courier', 9), anchor='w'
                          ).pack(anchor='w', padx=8, pady=1)
        if not any_warning:
            tk.Label(warn_frame, text='No major issues detected',
                      bg=biosim_theme('panel'), fg='#1a7f37',
                      font=('Arial', 10, 'italic')).pack(anchor='w', padx=8, pady=2)

        # Copy button
        tk.Button(card, text='Copy primer',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=lambda: (self.clipboard_clear(),
                                       self.clipboard_append(seq),
                                       self.status_lbl.config(
                                           text='Primer copied', fg='#1a7f37'))
                   ).pack(fill='x', padx=12, pady=(6, 10))

    def _do_primer_pair(self):
        """Analyze a manually entered primer pair (cross-dimer check)."""
        fw = ''.join(c.upper() for c in self.primer_fw_var.get() if c.isalpha())
        rv = ''.join(c.upper() for c in self.primer_rv_var.get() if c.isalpha())
        if not fw or not rv:
            self._primer_show_error('Enter both forward and reverse primers')
            return

        fw_result = score_primer(fw, tm_target=60)
        rv_result = score_primer(rv, tm_target=60)
        if not fw_result or not rv_result:
            self._primer_show_error('Primer length must be 15-35 nt')
            return

        cross = find_dimers(fw, rv, min_match=4)

        tm_diff = abs(fw_result['tm'] - rv_result['tm'])
        combined = (fw_result['score'] + rv_result['score']) / 2
        combined -= tm_diff * 2
        if cross and any(d['is_3_end'] for d in cross):
            combined -= 20

        # Build a fake pair object and render it
        pair = {
            'forward': fw_result,
            'reverse': rv_result,
            'product_length': '?',
            'tm_diff': tm_diff,
            'cross_dimers': cross,
            'combined_score': max(0, combined),
        }

        self._primer_clear_results()
        tk.Label(self.primer_results_inner,
                  text=('Анализ пары' if self.lang == 'ru' else 'Pair Analysis'),
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 14, 'bold')).pack(anchor='w', padx=12, pady=6)
        self._render_primer_pair_card(pair, 1)

    def _primer_show_error(self, msg):
        self._primer_clear_results()
        tk.Label(self.primer_results_inner, text='' + msg,
                  bg=biosim_theme('bg'), fg='#cf222e',
                  font=('Arial', 11, 'bold'),
                  padx=20, pady=20).pack(anchor='center')

    def _score_color(self, score):
        if score >= 85: return '#1a7f37'
        if score >= 70: return '#9a6700'
        if score >= 50: return '#bc4c00'
        return '#cf222e'

    def _tm_color(self, tm):
        if 55 <= tm <= 65: return '#1a7f37'
        if 50 <= tm <= 70: return '#9a6700'
        return '#cf222e'

    def _gc_color(self, gc):
        if 40 <= gc <= 60: return '#1a7f37'
        if 30 <= gc <= 70: return '#9a6700'
        return '#cf222e'

    def _clamp_color(self, status):
        return {'ok': '#1a7f37', 'acceptable': '#9a6700',
                'weak': '#bc4c00', 'too_high': '#cf222e',
                'too_short': '#cf222e'}.get(status, '#1c2128')

    def _build_tab_restriction(self):
        """Mol-3: Restriction sites finder, virtual digest, gel simulator, Gibson."""
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=self.T['tab_restriction'])
        self.tab_restriction = tab

        # ─── Top controls ───
        ctrl = tk.LabelFrame(tab,
                                text='  ' + ('Рестрикция и клонирование' if self.lang == 'ru'
                                                  else 'Restriction & Cloning') + '  ',
                                bg=biosim_theme('bg'), fg='#bc4c00',
                                font=('Arial', 11, 'bold'), relief='flat')
        ctrl.pack(fill='x', padx=4, pady=4)

        # Mode buttons
        mode_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        mode_frame.pack(fill='x', padx=8, pady=4)
        tk.Label(mode_frame, text=('Режим:' if self.lang == 'ru' else 'Mode:'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 10, 'bold')).pack(side='left', padx=(0, 8))

        self.restr_mode_var = tk.StringVar(value='map')
        for val, lbl, col in [
            ('map', '' + ('Карта сайтов' if self.lang == 'ru' else 'Restriction map'), '#0550ae'),
            ('digest', '' + ('Виртуальный digest' if self.lang == 'ru' else 'Virtual digest'), '#bc4c00'),
            ('gel', '' + ('Гель-электрофорез' if self.lang == 'ru' else 'Gel'), '#1a7f37'),
            ('gibson', 'Gibson assembly', '#8250df'),
        ]:
            tk.Radiobutton(mode_frame, text=lbl, value=val,
                            variable=self.restr_mode_var,
                            bg=biosim_theme('bg'), fg=col,
                            selectcolor='#eef1f5',
                            activebackground='#f5f7fa',
                            font=('Arial', 9, 'bold'),
                            command=self._switch_restr_mode
                            ).pack(side='left', padx=4)

        # Enzyme selector
        enz_row = tk.Frame(ctrl, bg=biosim_theme('bg'))
        enz_row.pack(fill='x', padx=8, pady=4)
        tk.Label(enz_row, text=('Ферменты:' if self.lang == 'ru' else 'Enzymes:'),
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 9, 'bold')).pack(side='left')

        # Selected enzymes (default: common ones)
        self.selected_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI']
        self.enz_display = tk.Label(enz_row,
                                       text=', '.join(self.selected_enzymes),
                                       bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                       font=('Courier', 9))
        self.enz_display.pack(side='left', padx=8)

        tk.Button(enz_row, text=('Выбрать...' if self.lang == 'ru' else 'Select...'),
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 9), cursor='hand2', padx=10, pady=2,
                   command=self._open_enzyme_selector
                   ).pack(side='left', padx=8)

        self.is_circular_var = tk.BooleanVar(value=False)
        tk.Checkbutton(enz_row, text=('Кольцевая (плазмида)' if self.lang == 'ru'
                                        else 'Circular (plasmid)'),
                        variable=self.is_circular_var,
                        bg=biosim_theme('bg'), fg='#9a6700',
                        selectcolor='#eef1f5',
                        activebackground='#f5f7fa',
                        font=('Arial', 9)).pack(side='left', padx=8)

        tk.Button(enz_row, text='▶ ' + ('Анализ' if self.lang == 'ru' else 'Analyze'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=3,
                   command=self._run_restriction
                   ).pack(side='right', padx=4)
        tk.Button(enz_row, text=('⭕ Карта плазмиды' if self.lang == 'ru'
                                  else '⭕ Plasmid map'),
                   bg=biosim_theme('panel'), highlightbackground='#8250df', fg='#8250df',
                   relief='solid', bd=1, font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=12, pady=3,
                   command=self._open_plasmid_map_dialog
                   ).pack(side='right', padx=4)

        # Gibson-specific input area (hidden by default)
        self.gibson_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        # Will be packed when gibson mode selected
        tk.Label(self.gibson_frame,
                  text=('Введи фрагменты (по одной строке, FASTA или DNA):'
                        if self.lang == 'ru' else
                        'Enter fragments (one per line, FASTA or DNA):'),
                  bg=biosim_theme('bg'), fg='#8250df',
                  font=('Arial', 9, 'bold')).pack(anchor='w', padx=4)
        gibson_inner = tk.Frame(self.gibson_frame, bg=biosim_theme('panel'))
        gibson_inner.pack(fill='x', padx=4, pady=2)
        g_sb = tk.Scrollbar(gibson_inner)
        g_sb.pack(side='right', fill='y')
        self.gibson_text = tk.Text(gibson_inner,
                                       bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                       insertbackground='#1c2128',
                                       font=('Courier', 9), relief='flat',
                                       height=4, padx=8, pady=4,
                                       yscrollcommand=g_sb.set)
        self.gibson_text.pack(side='left', fill='x', expand=True)
        g_sb.config(command=self.gibson_text.yview)

        # ─── Results area ───
        results_frame = tk.LabelFrame(tab,
                                          text='  ' + ('Результаты' if self.lang == 'ru'
                                                           else 'Results') + '  ',
                                          bg=biosim_theme('bg'), fg='#1a7f37',
                                          font=('Arial', 11, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=4, pady=4)

        canv = tk.Canvas(results_frame, bg=biosim_theme('bg'), highlightthickness=0)
        sb = tk.Scrollbar(results_frame, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        sb.pack(side='right', fill='y')
        self.restr_results_inner = tk.Frame(canv, bg=biosim_theme('bg'))
        self.restr_results_inner.bind('<Configure>',
                                          lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=self.restr_results_inner,
                             anchor='nw', width=1100)

        self._restr_show_initial_msg()

    def _switch_restr_mode(self):
        mode = self.restr_mode_var.get()
        if mode == 'gibson':
            self.gibson_frame.pack(fill='x', padx=8, pady=4)
        else:
            self.gibson_frame.pack_forget()
        self._restr_clear_results()
        self._restr_show_initial_msg()

    def _open_enzyme_selector(self):
        """Dialog to select restriction enzymes."""
        dlg = tk.Toplevel(self)
        dlg.title('Select Enzymes')
        dlg.geometry('540x600')
        dlg.configure(bg=biosim_theme('bg'))
        dlg.transient(self)

        tk.Label(dlg, text=('Выбери ферменты' if self.lang == 'ru'
                                       else 'Select restriction enzymes'),
                  bg=biosim_theme('bg'), fg='#bc4c00',
                  font=('Arial', 13, 'bold')).pack(pady=10)

        # Scrollable list with checkboxes
        canv = tk.Canvas(dlg, bg=biosim_theme('bg'), highlightthickness=0)
        sb = tk.Scrollbar(dlg, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True, padx=14, pady=4)
        sb.pack(side='right', fill='y')
        inner = tk.Frame(canv, bg=biosim_theme('bg'))
        inner.bind('<Configure>',
                    lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=inner, anchor='nw', width=500)

        enz_vars = {}
        for name in sorted(RESTRICTION_ENZYMES.keys()):
            enz = RESTRICTION_ENZYMES[name]
            var = tk.IntVar(value=1 if name in self.selected_enzymes else 0)
            enz_vars[name] = var
            row = tk.Frame(inner, bg=biosim_theme('panel'))
            row.pack(fill='x', pady=1, padx=2)
            tk.Checkbutton(row, text=name, variable=var,
                            bg=biosim_theme('panel'), fg='#0550ae',
                            selectcolor='#eef1f5',
                            activebackground='#ffffff',
                            font=('Arial', 10, 'bold'),
                            width=12, anchor='w').pack(side='left')
            tk.Label(row, text=enz['site'],
                      bg=biosim_theme('panel'), fg='#1a7f37',
                      font=('Courier', 10)).pack(side='left', padx=8)
            tk.Label(row, text='cut@' + str(enz['cut']),
                      bg=biosim_theme('panel'), fg=biosim_theme('text_dim'),
                      font=('Courier', 9)).pack(side='left', padx=4)
            tk.Label(row, text=enz['overhang'],
                      bg=biosim_theme('panel'), fg='#9a6700',
                      font=('Arial', 9, 'italic')).pack(side='left', padx=4)

        # Bottom buttons
        btn_row = tk.Frame(dlg, bg=biosim_theme('bg'))
        btn_row.pack(fill='x', padx=14, pady=10)

        def apply_sel():
            self.selected_enzymes = [n for n, v in enz_vars.items() if v.get()]
            self.enz_display.config(text=', '.join(self.selected_enzymes) or '(none)')
            dlg.destroy()

        def select_common():
            common = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'NotI', 'SalI']
            for n, v in enz_vars.items():
                v.set(1 if n in common else 0)

        def select_all():
            for v in enz_vars.values():
                v.set(1)

        def clear_all():
            for v in enz_vars.values():
                v.set(0)

        tk.Button(btn_row, text='Common', bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'),
                   relief='solid', bd=1, font=('Arial', 9), cursor='hand2',
                   padx=8, pady=3, command=select_common
                   ).pack(side='left', padx=2)
        tk.Button(btn_row, text='All', bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'),
                   relief='solid', bd=1, font=('Arial', 9), cursor='hand2',
                   padx=8, pady=3, command=select_all
                   ).pack(side='left', padx=2)
        tk.Button(btn_row, text='Clear', bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'),
                   relief='solid', bd=1, font=('Arial', 9), cursor='hand2',
                   padx=8, pady=3, command=clear_all
                   ).pack(side='left', padx=2)
        tk.Button(btn_row, text='Apply', bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=3, command=apply_sel
                   ).pack(side='right', padx=2)

    def _restr_clear_results(self):
        for w in self.restr_results_inner.winfo_children():
            w.destroy()

    def _restr_show_initial_msg(self):
        mode = self.restr_mode_var.get()
        msgs = {
            'map': ('Покажет все сайты выбранных ферментов в твоей ДНК'
                     if self.lang == 'ru' else
                     'Show all sites of selected enzymes in your DNA'),
            'digest': ('Разрежет ДНК ферментами и покажет фрагменты'
                        if self.lang == 'ru' else
                        'Cut DNA with enzymes and show fragments'),
            'gel': ('Симуляция гель-электрофореза для фрагментов'
                     if self.lang == 'ru' else
                     'Gel electrophoresis simulation'),
            'gibson': ('Стыковка фрагментов через homologous overlaps'
                         if self.lang == 'ru' else
                         'Stitch fragments via homologous overlaps'),
        }
        tk.Label(self.restr_results_inner, text=msgs.get(mode, ''),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 11, 'italic'),
                  padx=20, pady=20).pack(anchor='center')

    def _run_restriction(self):
        mode = self.restr_mode_var.get()
        if mode == 'map':
            self._do_restr_map()
        elif mode == 'digest':
            self._do_restr_digest()
        elif mode == 'gel':
            self._do_restr_gel()
        elif mode == 'gibson':
            self._do_gibson()

    def _do_restr_map(self):
        """Show all restriction sites in current sequence."""
        if not self.current_sequence:
            self._restr_show_error('Load a DNA sequence first')
            return
        if not self.selected_enzymes:
            self._restr_show_error('Select at least one enzyme')
            return

        sites = find_restriction_sites(self.current_sequence,
                                          enzymes=self.selected_enzymes)

        self._restr_clear_results()

        # Header
        tk.Label(self.restr_results_inner,
                  text=('Найдено сайтов' if self.lang == 'ru' else 'Sites found') +
                        ': {}'.format(len(sites)),
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 13, 'bold'), pady=4).pack(anchor='w', padx=12)

        if not sites:
            tk.Label(self.restr_results_inner,
                      text='No restriction sites found',
                      bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Arial', 11, 'italic')).pack(pady=10)
            return

        # Group by enzyme
        by_enzyme = {}
        for s in sites:
            by_enzyme.setdefault(s['enzyme'], []).append(s)

        # Summary table
        summary = tk.LabelFrame(self.restr_results_inner,
                                   text='  Summary  ',
                                   bg=biosim_theme('bg'), fg='#1a7f37',
                                   font=('Arial', 10, 'bold'), relief='flat')
        summary.pack(fill='x', padx=10, pady=4)

        for enz_name in sorted(by_enzyme.keys()):
            sites_for_enz = by_enzyme[enz_name]
            row = tk.Frame(summary, bg=biosim_theme('bg'))
            row.pack(fill='x', padx=8, pady=2)
            enz = RESTRICTION_ENZYMES[enz_name]
            tk.Label(row, text=enz_name, bg=biosim_theme('bg'), fg='#0550ae',
                      font=('Arial', 10, 'bold'), width=10, anchor='w'
                      ).pack(side='left')
            tk.Label(row, text=enz['site'], bg=biosim_theme('bg'), fg='#1a7f37',
                      font=('Courier', 10), width=12, anchor='w'
                      ).pack(side='left')
            tk.Label(row, text='{} sites'.format(len(sites_for_enz)),
                      bg=biosim_theme('bg'), fg='#9a6700',
                      font=('Courier', 10), width=12, anchor='w'
                      ).pack(side='left')
            positions = ', '.join(str(s['position'] + 1) for s in sites_for_enz[:8])
            if len(sites_for_enz) > 8:
                positions += '...'
            tk.Label(row, text='@ ' + positions,
                      bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Courier', 9), anchor='w'
                      ).pack(side='left', padx=4)

        # Visual map
        self._draw_restriction_map(by_enzyme)

    def _draw_restriction_map(self, sites_by_enz):
        """Draw a visual linear restriction map."""
        map_frame = tk.LabelFrame(self.restr_results_inner,
                                     text='  Visual map  ',
                                     bg=biosim_theme('bg'), fg='#0969da',
                                     font=('Arial', 10, 'bold'), relief='flat')
        map_frame.pack(fill='x', padx=10, pady=4)

        canv_w = 1000
        canv_h = max(120, 50 + len(sites_by_enz) * 30)
        c = tk.Canvas(map_frame, bg=biosim_theme('panel'), height=canv_h, width=canv_w,
                        highlightthickness=0)
        c.pack(fill='x', padx=4, pady=4)

        seq_len = len(self.current_sequence)
        margin = 30
        plot_w = canv_w - 2 * margin

        # Main sequence line
        center_y = canv_h // 2
        c.create_rectangle(margin, center_y - 4, margin + plot_w, center_y + 4,
                              fill='#d0d7de', outline='')

        # Tick marks
        n_ticks = 10
        for i in range(n_ticks + 1):
            x = margin + plot_w * i / n_ticks
            pos = int(seq_len * i / n_ticks)
            c.create_line(x, center_y - 4, x, center_y + 8, fill='#57606a')
            c.create_text(x, center_y + 16,
                            text=str(pos), fill='#57606a',
                            font=('Courier', 8))

        # Plot each enzyme on its own track above/below
        colors = ['#1a7f37', '#0550ae', '#9a6700', '#bc4c00',
                   '#8250df', '#cf222e', '#8250df', '#bc4c00']
        y_offset = 30
        for i, (enz_name, sites) in enumerate(sorted(sites_by_enz.items())):
            color = colors[i % len(colors)]
            y_label = center_y - y_offset - (i % 4) * 16
            c.create_text(margin - 5, y_label, text=enz_name,
                            fill=color, font=('Arial', 8, 'bold'),
                            anchor='e')
            for s in sites:
                x = margin + plot_w * s['position'] / max(seq_len, 1)
                c.create_line(x, center_y, x, y_label + 8,
                                fill=color, width=2)
                c.create_oval(x - 3, y_label + 5, x + 3, y_label + 11,
                                fill=color, outline='')

    def _do_restr_digest(self):
        """Virtual digest — show resulting fragments."""
        if not self.current_sequence:
            self._restr_show_error('Load a DNA sequence first')
            return
        if not self.selected_enzymes:
            self._restr_show_error('Select at least one enzyme')
            return

        fragments = digest_sequence(self.current_sequence,
                                       enzymes=self.selected_enzymes,
                                       is_circular=self.is_circular_var.get())

        self._restr_clear_results()

        # Header
        circular_text = '(circular)' if self.is_circular_var.get() else '(linear)'
        tk.Label(self.restr_results_inner,
                  text=('Получено фрагментов' if self.lang == 'ru'
                                  else 'Fragments produced') +
                        ': {} {}'.format(len(fragments), circular_text),
                  bg=biosim_theme('bg'), fg='#bc4c00',
                  font=('Arial', 13, 'bold'), pady=4).pack(anchor='w', padx=12)

        # List of fragments
        for i, frag in enumerate(fragments, 1):
            card = tk.Frame(self.restr_results_inner, bg=biosim_theme('panel'),
                              highlightthickness=2,
                              highlightbackground=biosim_theme('border'))
            card.pack(fill='x', padx=10, pady=4)

            head = tk.Frame(card, bg=biosim_theme('panel'))
            head.pack(fill='x', padx=10, pady=(8, 2))
            tk.Label(head, text='Fragment #{}'.format(i),
                      bg=biosim_theme('panel'), fg='#0550ae',
                      font=('Arial', 11, 'bold')).pack(side='left')
            tk.Label(head, text='  {} bp'.format(frag['length']),
                      bg=biosim_theme('panel'), fg='#1a7f37',
                      font=('Courier', 11, 'bold')).pack(side='left', padx=8)

            le = frag.get('left_enzyme')
            re = frag.get('right_enzyme')
            ends_text = '['
            ends_text += le if le else 'start'
            ends_text += ' — '
            ends_text += re if re else 'end'
            ends_text += ']'
            tk.Label(head, text=ends_text,
                      bg=biosim_theme('panel'), fg='#9a6700',
                      font=('Courier', 9)).pack(side='left', padx=8)

            # Sequence preview (first 60 + last 30 bp)
            seq = frag['sequence']
            if len(seq) > 100:
                preview = seq[:60] + '   ...   ' + seq[-30:]
            else:
                preview = seq
            tk.Label(card, text=preview,
                      bg=biosim_theme('panel'), fg=biosim_theme('text'),
                      font=('Courier', 9), wraplength=1000,
                      justify='left').pack(anchor='w', padx=14, pady=(0, 8))

            def make_copy(s):
                def copy():
                    self.clipboard_clear()
                    self.clipboard_append(s)
                    self.status_lbl.config(text='Fragment copied', fg='#1a7f37')
                return copy

            tk.Button(card, text='Copy sequence',
                       bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                       font=('Arial', 9, 'bold'), cursor='hand2',
                       padx=10, pady=2,
                       command=make_copy(seq)
                       ).pack(anchor='w', padx=14, pady=(0, 8))

    def _do_restr_gel(self):
        """Simulate gel electrophoresis."""
        if not self.current_sequence:
            self._restr_show_error('Load a DNA sequence first')
            return
        if not self.selected_enzymes:
            self._restr_show_error('Select at least one enzyme')
            return

        fragments = digest_sequence(self.current_sequence,
                                       enzymes=self.selected_enzymes,
                                       is_circular=self.is_circular_var.get())

        self._restr_clear_results()

        tk.Label(self.restr_results_inner,
                  text=('Гель электрофорез' if self.lang == 'ru'
                                  else 'Gel electrophoresis'),
                  bg=biosim_theme('bg'), fg='#1a7f37',
                  font=('Arial', 13, 'bold'), pady=4).pack(anchor='w', padx=12)

        # Gel canvas
        gel_w = 600
        gel_h = 500
        c = tk.Canvas(self.restr_results_inner, bg=biosim_theme('panel'),
                        width=gel_w, height=gel_h,
                        highlightthickness=2,
                        highlightbackground=biosim_theme('border'))
        c.pack(padx=10, pady=8)

        # Background gradient (gel appearance)
        for y in range(gel_h):
            shade = int(20 + (y / gel_h) * 5)
            color = '#{:02x}{:02x}{:02x}'.format(shade + 8, shade + 8, shade + 12)
            c.create_line(0, y, gel_w, y, fill=color)

        # Wells (top)
        well_w = 60
        n_lanes = 3  # ladder, sample, marker
        for lane in range(n_lanes):
            x = 80 + lane * (well_w + 60)
            c.create_rectangle(x, 10, x + well_w, 30,
                                  fill='#eef1f5', outline='#f5f7fa', width=2)
            labels = ['Ladder', 'Your DNA', '']
            c.create_text(x + well_w / 2, 5, text=labels[lane],
                            fill='#1c2128', font=('Arial', 9, 'bold'),
                            anchor='s')

        # Ladder bands (lane 0)
        ladder_x = 80
        ladder = DNA_LADDERS['1kb_plus']
        for bp in ladder:
            y = 30 + estimate_gel_migration(bp, gel_pct=1.0) * (gel_h - 60)
            c.create_rectangle(ladder_x + 5, y - 1, ladder_x + well_w - 5, y + 2,
                                  fill='#1c2128', outline='')
            # Label on the right
            c.create_text(ladder_x + well_w + 5, y,
                            text='{}'.format(bp) if bp >= 1000 else '{}'.format(bp),
                            fill='#57606a', font=('Courier', 8),
                            anchor='w')

        # Sample bands (lane 1)
        sample_x = 80 + (well_w + 60)
        for frag in fragments:
            bp = frag['length']
            y = 30 + estimate_gel_migration(bp, gel_pct=1.0) * (gel_h - 60)
            # Brightness proportional to amount (assume equimolar — same brightness)
            c.create_rectangle(sample_x + 5, y - 2, sample_x + well_w - 5, y + 3,
                                  fill='#1a7f37', outline='')
            # Label
            c.create_text(sample_x + well_w + 5, y,
                            text='{} bp'.format(bp),
                            fill='#1a7f37', font=('Courier', 8, 'bold'),
                            anchor='w')

        # Direction arrow
        c.create_text(20, gel_h // 2, text='—\n\n+',
                        fill='#57606a', font=('Arial', 10, 'bold'),
                        anchor='w')

        # Info
        tk.Label(self.restr_results_inner,
                  text='Lane 1: 1kb+ ladder · Lane 2: digest of your DNA',
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9, 'italic')).pack(anchor='w', padx=14)
        tk.Label(self.restr_results_inner,
                  text=('1% agarose, 1× TAE buffer, ~80V, 30 min (simulation)'
                        if self.lang == 'en' else
                        '1% агароза, 1× TAE, ~80V, 30 мин (симуляция)'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9, 'italic')).pack(anchor='w', padx=14)

    def _do_gibson(self):
        """Gibson assembly of fragments."""
        text = self.gibson_text.get('1.0', 'end').strip()
        if not text:
            self._restr_show_error('Enter fragments in the box above')
            return

        # Parse fragments — one per line, allow FASTA
        fragments = []
        current = []
        for line in text.split('\n'):
            line = line.strip()
            if not line:
                if current:
                    fragments.append(''.join(current).upper())
                    current = []
                continue
            if line.startswith('>'):
                if current:
                    fragments.append(''.join(current).upper())
                    current = []
                continue
            # Clean — only DNA letters
            cleaned = ''.join(c for c in line.upper() if c in 'ATGCN')
            current.append(cleaned)
        if current:
            fragments.append(''.join(current))

        fragments = [f for f in fragments if len(f) >= 20]

        if len(fragments) < 2:
            self._restr_show_error('Need at least 2 fragments of ≥20 bp each')
            return

        result = gibson_assembly(fragments, overlap_min=15, overlap_max=40)

        self._restr_clear_results()
        if not result['success']:
            tk.Label(self.restr_results_inner,
                      text='' + result['message'],
                      bg=biosim_theme('bg'), fg='#cf222e',
                      font=('Arial', 12, 'bold'),
                      pady=10).pack(padx=14)
            # Show what overlaps were found
            if result['overlaps']:
                tk.Label(self.restr_results_inner,
                          text='Overlaps found before failure:',
                          bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                          font=('Arial', 10, 'italic')).pack(anchor='w', padx=14)
                for ov in result['overlaps']:
                    tk.Label(self.restr_results_inner,
                              text='  Fragment {} {}: {} bp overlap'.format(
                                  ov['frag_i'], ov['frag_j'], ov['length']),
                              bg=biosim_theme('bg'), fg='#0550ae',
                              font=('Courier', 9)).pack(anchor='w', padx=20)
            return

        # Success
        tk.Label(self.restr_results_inner,
                  text='' + result['message'],
                  bg=biosim_theme('bg'), fg='#1a7f37',
                  font=('Arial', 13, 'bold'), pady=6).pack(anchor='w', padx=12)

        # Overlaps table
        ov_frame = tk.LabelFrame(self.restr_results_inner,
                                    text='  Overlaps  ',
                                    bg=biosim_theme('bg'), fg='#8250df',
                                    font=('Arial', 10, 'bold'), relief='flat')
        ov_frame.pack(fill='x', padx=10, pady=4)
        for ov in result['overlaps']:
            row = tk.Frame(ov_frame, bg=biosim_theme('bg'))
            row.pack(fill='x', padx=8, pady=2)
            tk.Label(row,
                      text='Fragment {} {}'.format(ov['frag_i'], ov['frag_j']),
                      bg=biosim_theme('bg'), fg='#0550ae',
                      font=('Arial', 10, 'bold'), width=22, anchor='w'
                      ).pack(side='left')
            tk.Label(row, text='{} bp:'.format(ov['length']),
                      bg=biosim_theme('bg'), fg='#9a6700',
                      font=('Courier', 10)).pack(side='left', padx=4)
            tk.Label(row, text=ov['overlap_seq'][:50],
                      bg=biosim_theme('bg'), fg='#1a7f37',
                      font=('Courier', 9)).pack(side='left', padx=4)

        # Assembled sequence
        assembled = result['assembled']
        result_frame = tk.LabelFrame(self.restr_results_inner,
                                        text='  Assembled product · {} bp  '.format(
                                            len(assembled)),
                                        bg=biosim_theme('bg'), fg='#1a7f37',
                                        font=('Arial', 10, 'bold'), relief='flat')
        result_frame.pack(fill='both', expand=True, padx=10, pady=4)

        result_text = tk.Text(result_frame, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                  font=('Courier', 9), relief='flat',
                                  wrap='word', height=8, padx=10, pady=8)
        result_text.pack(fill='both', expand=True, padx=4, pady=4)
        result_text.insert('1.0', format_sequence_blocks(assembled,
                                                              block_size=10,
                                                              line_size=60))
        result_text.config(state='disabled')

        # Copy + Load buttons
        btn_row = tk.Frame(result_frame, bg=biosim_theme('bg'))
        btn_row.pack(fill='x', padx=4, pady=4)

        def copy_assembled():
            self.clipboard_clear()
            self.clipboard_append(assembled)
            self.status_lbl.config(text='Assembled product copied',
                                       fg='#1a7f37')

        def load_assembled():
            self.current_sequence = assembled
            self.sequence_name = 'Gibson_product'
            self.sequence_type = 'dna'
            if hasattr(self, 'seq_text'):
                self.seq_text.delete('1.0', 'end')
                self.seq_text.insert('1.0', assembled)
                if hasattr(self, 'seq_name_var'):
                    self.seq_name_var.set('Gibson_product')
                self._refresh_seq_analysis()
                self._refresh_seq_highlights()
                self.notebook.select(self.tab_sequence)
                self.status_lbl.config(
                    text='Loaded {} bp into Sequence Editor'.format(
                        len(assembled)),
                    fg='#1a7f37')

        tk.Button(btn_row, text='Copy', bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'),
                   relief='solid', bd=1, font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=12, pady=3, command=copy_assembled
                   ).pack(side='left', padx=2)
        tk.Button(btn_row, text='Load into editor',
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=12, pady=3, command=load_assembled
                   ).pack(side='left', padx=2)

    def _restr_show_error(self, msg):
        self._restr_clear_results()
        tk.Label(self.restr_results_inner, text='' + msg,
                  bg=biosim_theme('bg'), fg='#cf222e',
                  font=('Arial', 11, 'bold'),
                  padx=20, pady=20).pack(anchor='center')

    def _build_tab_alignment(self):
        """Mol-4: Pairwise alignment + BLAST-like search."""
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=self.T['tab_align'])
        self.tab_alignment = tab

        # ─── Top controls ───
        ctrl = tk.LabelFrame(tab,
                                text='  ' + ('Выравнивание' if self.lang == 'ru'
                                                  else 'Sequence Alignment') + '  ',
                                bg=biosim_theme('bg'), fg='#0969da',
                                font=('Arial', 11, 'bold'), relief='flat')
        ctrl.pack(fill='x', padx=4, pady=4)

        # Mode buttons
        mode_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        mode_frame.pack(fill='x', padx=8, pady=4)
        tk.Label(mode_frame, text=('Режим:' if self.lang == 'ru' else 'Mode:'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 10, 'bold')).pack(side='left', padx=(0, 8))

        self.align_mode_var = tk.StringVar(value='global')
        for val, lbl, col in [
            ('global', '' + ('Глобальное (Needleman-Wunsch)' if self.lang == 'ru' else 'Global (Needleman-Wunsch)'), '#1a7f37'),
            ('local',  '' + ('Локальное (Smith-Waterman)' if self.lang == 'ru' else 'Local (Smith-Waterman)'), '#0969da'),
            ('blast',  'BLAST-like vs built-in DB', '#8250df'),
        ]:
            tk.Radiobutton(mode_frame, text=lbl, value=val,
                            variable=self.align_mode_var,
                            bg=biosim_theme('bg'), fg=col,
                            selectcolor='#eef1f5',
                            activebackground='#f5f7fa',
                            font=('Arial', 9, 'bold'),
                            command=self._switch_align_mode
                            ).pack(side='left', padx=4)

        # Sequence inputs
        seqs_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        seqs_frame.pack(fill='x', padx=8, pady=4)

        tk.Label(seqs_frame, text=('Последовательность 1:' if self.lang == 'ru'
                                       else 'Sequence 1:'),
                  bg=biosim_theme('bg'), fg='#1a7f37',
                  font=('Arial', 9, 'bold')).pack(anchor='w')
        s1_inner = tk.Frame(seqs_frame, bg=biosim_theme('panel'))
        s1_inner.pack(fill='x', pady=2)
        self.align_seq1_text = tk.Text(s1_inner, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                            insertbackground='#1c2128',
                                            font=('Courier', 9), relief='flat',
                                            height=3, padx=8, pady=4,
                                            wrap='word')
        self.align_seq1_text.pack(fill='x', expand=True)
        # Use current sequence button
        tk.Button(s1_inner, text='Use current',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 8, 'bold'), cursor='hand2',
                   padx=6, pady=1,
                   command=self._align_use_current_for_1
                   ).pack(side='right', padx=2, pady=2)

        # Sequence 2 (hidden in BLAST mode)
        self.align_seq2_frame = tk.Frame(seqs_frame, bg=biosim_theme('bg'))
        tk.Label(self.align_seq2_frame, text=('Последовательность 2:' if self.lang == 'ru'
                                                  else 'Sequence 2:'),
                  bg=biosim_theme('bg'), fg='#bc4c00',
                  font=('Arial', 9, 'bold')).pack(anchor='w')
        s2_inner = tk.Frame(self.align_seq2_frame, bg=biosim_theme('panel'))
        s2_inner.pack(fill='x', pady=2)
        self.align_seq2_text = tk.Text(s2_inner, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                            insertbackground='#1c2128',
                                            font=('Courier', 9), relief='flat',
                                            height=3, padx=8, pady=4,
                                            wrap='word')
        self.align_seq2_text.pack(fill='x', expand=True)
        self.align_seq2_frame.pack(fill='x', pady=2)

        # Parameters row
        param_row = tk.Frame(ctrl, bg=biosim_theme('bg'))
        param_row.pack(fill='x', padx=8, pady=4)

        self.align_protein_var = tk.BooleanVar(value=False)
        tk.Checkbutton(param_row, text=('Белок (BLOSUM62)' if self.lang == 'ru'
                                          else 'Protein (BLOSUM62)'),
                        variable=self.align_protein_var,
                        bg=biosim_theme('bg'), fg='#8250df',
                        selectcolor='#eef1f5',
                        activebackground='#f5f7fa',
                        font=('Arial', 9, 'bold')).pack(side='left', padx=4)

        tk.Label(param_row, text='Match:', bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9)).pack(side='left', padx=(10, 2))
        self.align_match_var = tk.IntVar(value=2)
        tk.Entry(param_row, textvariable=self.align_match_var,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=4
                  ).pack(side='left', padx=2)

        tk.Label(param_row, text='Mismatch:', bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9)).pack(side='left', padx=(8, 2))
        self.align_mismatch_var = tk.IntVar(value=-1)
        tk.Entry(param_row, textvariable=self.align_mismatch_var,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=4
                  ).pack(side='left', padx=2)

        tk.Label(param_row, text='Gap:', bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9)).pack(side='left', padx=(8, 2))
        self.align_gap_var = tk.IntVar(value=-2)
        tk.Entry(param_row, textvariable=self.align_gap_var,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Arial', 9), width=4
                  ).pack(side='left', padx=2)

        tk.Button(param_row, text='▶ ' + ('Выровнять' if self.lang == 'ru' else 'Align'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=16, pady=4,
                   command=self._run_alignment
                   ).pack(side='right', padx=4)
        tk.Button(param_row, text=('🔍 Варианты (SNP)' if self.lang == 'ru'
                                    else '🔍 Variants (SNP)'),
                   bg=biosim_theme('panel'), highlightbackground='#cf222e', fg='#cf222e',
                   relief='solid', bd=1, font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=12, pady=4,
                   command=self._open_variants_dialog
                   ).pack(side='right', padx=4)

        # ─── Results area ───
        results_frame = tk.LabelFrame(tab,
                                          text='  ' + ('Результаты' if self.lang == 'ru'
                                                           else 'Results') + '  ',
                                          bg=biosim_theme('bg'), fg='#1a7f37',
                                          font=('Arial', 11, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=4, pady=4)

        canv = tk.Canvas(results_frame, bg=biosim_theme('bg'), highlightthickness=0)
        sb = tk.Scrollbar(results_frame, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        sb.pack(side='right', fill='y')
        self.align_results_inner = tk.Frame(canv, bg=biosim_theme('bg'))
        self.align_results_inner.bind('<Configure>',
                                            lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=self.align_results_inner,
                             anchor='nw', width=1100)

        self._align_show_initial_msg()

    def _switch_align_mode(self):
        mode = self.align_mode_var.get()
        if mode == 'blast':
            self.align_seq2_frame.pack_forget()
        else:
            self.align_seq2_frame.pack(fill='x', pady=2)
        self._align_clear_results()
        self._align_show_initial_msg()

    def _align_use_current_for_1(self):
        """Load current_sequence into sequence 1 input."""
        if self.current_sequence:
            self.align_seq1_text.delete('1.0', 'end')
            self.align_seq1_text.insert('1.0', self.current_sequence)
            if self.sequence_type == 'protein':
                self.align_protein_var.set(True)

    def _align_clear_results(self):
        for w in self.align_results_inner.winfo_children():
            w.destroy()

    def _align_show_initial_msg(self):
        mode = self.align_mode_var.get()
        msgs = {
            'global': ('Выровняй две последовательности полностью (от начала до конца)'
                        if self.lang == 'ru' else
                        'Align two sequences fully (end-to-end)'),
            'local':  ('Найди лучший локальный участок выравнивания'
                        if self.lang == 'ru' else
                        'Find best local alignment region'),
            'blast':  ('Поиск похожих в библиотеке белков (GFP, p53, LacI, и т.д.)'
                        if self.lang == 'ru' else
                        'Search built-in protein library (GFP, p53, LacI, etc.)'),
        }
        tk.Label(self.align_results_inner, text=msgs.get(mode, ''),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 11, 'italic'),
                  padx=20, pady=20).pack(anchor='center')

    def _run_alignment(self):
        mode = self.align_mode_var.get()
        # Get sequence 1
        seq1 = ''.join(c for c in self.align_seq1_text.get('1.0', 'end').upper()
                        if c.isalpha())
        if not seq1:
            self._align_show_error('Enter Sequence 1' if self.lang == 'en'
                                       else 'Введи Последовательность 1')
            return

        if mode == 'blast':
            self._do_blast(seq1)
        else:
            seq2 = ''.join(c for c in self.align_seq2_text.get('1.0', 'end').upper()
                            if c.isalpha())
            if not seq2:
                self._align_show_error('Enter Sequence 2')
                return
            self._do_pairwise_align(seq1, seq2, local=(mode == 'local'))

    def _do_pairwise_align(self, seq1, seq2, local=False):
        """Run pairwise alignment and display result."""
        is_prot = self.align_protein_var.get()
        match = self.align_match_var.get()
        mismatch = self.align_mismatch_var.get()
        gap = self.align_gap_var.get()

        self.status_lbl.config(text=('Aligning...' if self.lang == 'en'
                                          else 'Выравниваю...'),
                                    fg='#0969da')
        self.update()

        try:
            if local:
                result = smith_waterman(seq1, seq2, match=match, mismatch=mismatch,
                                            gap=gap, is_protein=is_prot)
            else:
                result = needleman_wunsch(seq1, seq2, match=match, mismatch=mismatch,
                                              gap=gap, is_protein=is_prot)
        except Exception as e:
            self._align_show_error('Error: ' + str(e)[:80])
            biosim_log('error', 'alignment: ' + str(e))
            return

        if result is None:
            self._align_show_error('Alignment failed (sequence too long? max 5000 nt)')
            return

        self._align_clear_results()
        self._render_alignment_result(result, local=local,
                                          name1='Sequence 1', name2='Sequence 2')
        self.status_lbl.config(
            text='' + 'Aligned' + ' · Identity: {:.1f}%'.format(result['identity']),
            fg='#1a7f37')

    def _render_alignment_result(self, result, local=False, name1='Seq1', name2='Seq2'):
        """Display alignment with monospace pairing view."""
        # Header card with stats
        head = tk.LabelFrame(self.align_results_inner,
                                text='  ' + ('Статистика' if self.lang == 'ru'
                                                  else 'Statistics') + '  ',
                                bg=biosim_theme('bg'), fg='#1a7f37',
                                font=('Arial', 10, 'bold'), relief='flat')
        head.pack(fill='x', padx=10, pady=4)

        stats = tk.Frame(head, bg=biosim_theme('bg'))
        stats.pack(fill='x', padx=10, pady=6)

        items = [
            ('Type:', 'Local' if local else 'Global'),
            ('Score:', '{}'.format(result['score'])),
            ('Length:', '{}'.format(result['length'])),
            ('Identity:', '{:.1f}% ({}/{})'.format(
                result['identity'], result['matches'], result['length'])),
            ('Similarity:', '{:.1f}%'.format(result['similarity'])),
            ('Gaps:', '{}'.format(result['gaps'])),
        ]
        if local:
            items.append(('Region 1:', '{}-{}'.format(result.get('start1', 0) + 1,
                                                          result.get('end1', 0))))
            items.append(('Region 2:', '{}-{}'.format(result.get('start2', 0) + 1,
                                                          result.get('end2', 0))))
        for i, (k, v) in enumerate(items):
            col = i // 3
            row = i % 3
            tk.Label(stats, text=k, bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Arial', 9), anchor='e', width=12
                      ).grid(row=row, column=col * 2, sticky='e', padx=4, pady=1)
            color = '#1a7f37' if 'Identity' in k or 'Similarity' in k else '#1c2128'
            tk.Label(stats, text=v, bg=biosim_theme('bg'), fg=color,
                      font=('Courier', 10, 'bold'), anchor='w', width=20
                      ).grid(row=row, column=col * 2 + 1, sticky='w', padx=4, pady=1)

        # Alignment view (monospace 3-line format)
        aln_frame = tk.LabelFrame(self.align_results_inner,
                                      text='  ' + ('Выравнивание' if self.lang == 'ru'
                                                       else 'Alignment') + '  ',
                                      bg=biosim_theme('bg'), fg='#0550ae',
                                      font=('Arial', 10, 'bold'), relief='flat')
        aln_frame.pack(fill='both', expand=True, padx=10, pady=4)

        aln_inner = tk.Frame(aln_frame, bg=biosim_theme('panel'))
        aln_inner.pack(fill='both', expand=True, padx=4, pady=4)
        aln_sb = tk.Scrollbar(aln_inner)
        aln_sb.pack(side='right', fill='y')
        aln_txt = tk.Text(aln_inner, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                              font=('Courier', 10), relief='flat',
                              wrap='none', height=14, padx=10, pady=8,
                              yscrollcommand=aln_sb.set)
        aln_txt.pack(side='left', fill='both', expand=True)
        aln_sb.config(command=aln_txt.yview)

        # Format: line1 = seq1, line2 = match indicators, line3 = seq2
        # Break into chunks of 60
        a1 = result['aligned1']
        a2 = result['aligned2']
        is_prot = self.align_protein_var.get()
        chunk_size = 60

        # Tags for highlighting
        aln_txt.tag_config('match', foreground='#1a7f37',
                              font=('Courier', 10, 'bold'))
        aln_txt.tag_config('similar', foreground='#9a6700',
                              font=('Courier', 10))
        aln_txt.tag_config('mismatch', foreground='#cf222e',
                              font=('Courier', 10))
        aln_txt.tag_config('gap', foreground=biosim_theme('text_dim'),
                              font=('Courier', 10))
        aln_txt.tag_config('label', foreground='#0550ae',
                              font=('Courier', 10, 'bold'))
        aln_txt.tag_config('pos', foreground=biosim_theme('text_dim'),
                              font=('Courier', 9))

        # Compute starting positions (for local alignment, use start1/start2; else 0)
        pos1 = result.get('start1', 0)
        pos2 = result.get('start2', 0)

        for chunk_start in range(0, len(a1), chunk_size):
            chunk1 = a1[chunk_start:chunk_start + chunk_size]
            chunk2 = a2[chunk_start:chunk_start + chunk_size]

            # Count actual sequence chars (excluding gaps) for position labels
            real1_start = pos1 + 1
            real2_start = pos2 + 1
            
            # Line 1: Seq1
            aln_txt.insert('end', '{} {:>5} '.format(name1[:6], real1_start), 'label')
            for c in chunk1:
                tag = 'gap' if c == '-' else ''
                aln_txt.insert('end', c, tag)
            # Count chars consumed
            pos1 += sum(1 for c in chunk1 if c != '-')
            aln_txt.insert('end', ' {:>5}\n'.format(pos1), 'pos')

            # Line 2: match indicators
            aln_txt.insert('end', '       ' + ' ' * 6)
            for c1, c2 in zip(chunk1, chunk2):
                if c1 == '-' or c2 == '-':
                    aln_txt.insert('end', ' ', 'gap')
                elif c1 == c2:
                    aln_txt.insert('end', '|', 'match')
                elif is_prot and blosum62_score(c1, c2) > 0:
                    aln_txt.insert('end', ':', 'similar')
                else:
                    aln_txt.insert('end', '.', 'mismatch')
            aln_txt.insert('end', '\n')

            # Line 3: Seq2
            aln_txt.insert('end', '{} {:>5} '.format(name2[:6], real2_start), 'label')
            for c in chunk2:
                tag = 'gap' if c == '-' else ''
                aln_txt.insert('end', c, tag)
            pos2 += sum(1 for c in chunk2 if c != '-')
            aln_txt.insert('end', ' {:>5}\n\n'.format(pos2), 'pos')

        aln_txt.config(state='disabled')

        # Copy button
        def copy_alignment():
            text = '>{}\n{}\n>{}\n{}'.format(name1, a1, name2, a2)
            self.clipboard_clear()
            self.clipboard_append(text)
            self.status_lbl.config(text='Alignment copied', fg='#1a7f37')

        tk.Button(aln_frame, text='Copy alignment (FASTA)',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=10, pady=3, command=copy_alignment
                   ).pack(anchor='e', padx=8, pady=4)

    def _do_blast(self, query):
        """BLAST-like search against built-in protein DB."""
        is_prot = self.align_protein_var.get()
        if not is_prot:
            # Auto-enable for built-in DB
            self.align_protein_var.set(True)
            is_prot = True

        self.status_lbl.config(text=('Searching DB...' if self.lang == 'en'
                                          else 'Поиск в БД...'),
                                    fg='#0969da')
        self.update()

        try:
            results = blast_like_search(query, BUILT_IN_PROTEIN_DB,
                                            word_size=4,  # shorter words for protein
                                            max_results=10,
                                            is_protein=True)
        except Exception as e:
            self._align_show_error('BLAST error: ' + str(e)[:80])
            biosim_log('error', 'blast: ' + str(e))
            return

        self._align_clear_results()

        if not results:
            tk.Label(self.align_results_inner,
                      text=('Совпадений не найдено' if self.lang == 'ru'
                                       else 'No hits found'),
                      bg=biosim_theme('bg'), fg='#cf222e',
                      font=('Arial', 12, 'bold'),
                      pady=20).pack()
            return

        # Header
        tk.Label(self.align_results_inner,
                  text=('Найдено совпадений' if self.lang == 'ru'
                                  else 'Hits found') + ': {}'.format(len(results)),
                  bg=biosim_theme('bg'), fg='#8250df',
                  font=('Arial', 13, 'bold'), pady=4).pack(anchor='w', padx=12)

        # Each hit as a card
        for i, hit in enumerate(results, 1):
            card = tk.Frame(self.align_results_inner, bg=biosim_theme('panel'),
                              highlightthickness=2,
                              highlightbackground=self._blast_score_color(hit['identity']))
            card.pack(fill='x', padx=10, pady=4)

            head = tk.Frame(card, bg=biosim_theme('panel'))
            head.pack(fill='x', padx=10, pady=(8, 2))
            tk.Label(head, text='#{} {}'.format(i, hit['name']),
                      bg=biosim_theme('panel'), fg='#0550ae',
                      font=('Arial', 12, 'bold')).pack(side='left')
            tk.Label(head,
                      text='  Score: {} · Identity: {:.1f}% · E-value: {:.2e}'.format(
                          hit['score'], hit['identity'], hit['e_value']),
                      bg=biosim_theme('panel'), fg='#9a6700',
                      font=('Courier', 9)).pack(side='left', padx=8)

            if hit.get('description'):
                tk.Label(card, text=hit['description'],
                          bg=biosim_theme('panel'), fg=biosim_theme('text_dim'),
                          font=('Arial', 9, 'italic'),
                          wraplength=900, justify='left',
                          ).pack(anchor='w', padx=14, pady=2)

            # Show alignment preview (just first chunk)
            aln = hit['alignment']
            preview1 = aln['aligned1'][:60]
            preview2 = aln['aligned2'][:60]
            preview_text = preview1 + '\n' + preview2
            tk.Label(card, text=preview_text,
                      bg=biosim_theme('panel'), fg=biosim_theme('text'),
                      font=('Courier', 9), justify='left',
                      anchor='w').pack(anchor='w', padx=14, pady=2)

            # View full alignment button
            def make_view(h):
                def view_full():
                    self._show_blast_full_alignment(h, query)
                return view_full

            tk.Button(card, text='View full alignment',
                       bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                       font=('Arial', 9, 'bold'), cursor='hand2',
                       padx=10, pady=2,
                       command=make_view(hit)
                       ).pack(anchor='w', padx=14, pady=(2, 8))

    def _show_blast_full_alignment(self, hit, query):
        """Show full alignment of BLAST hit in a separate dialog."""
        dlg = tk.Toplevel(self)
        dlg.title('Full alignment: ' + hit['name'])
        dlg.geometry('900x600')
        dlg.configure(bg=biosim_theme('bg'))
        dlg.transient(self)

        tk.Label(dlg, text='' + hit['name'],
                  bg=biosim_theme('bg'), fg='#8250df',
                  font=('Arial', 14, 'bold')).pack(pady=8)
        tk.Label(dlg, text=hit.get('description', ''),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 10, 'italic')).pack()
        tk.Label(dlg,
                  text='Identity: {:.1f}% · Length: {} · E-value: {:.2e}'.format(
                      hit['identity'], hit['length'], hit['e_value']),
                  bg=biosim_theme('bg'), fg='#1a7f37',
                  font=('Courier', 10)).pack(pady=4)

        # Render alignment
        # Save inner & make inner just for this dialog
        outer = tk.Frame(dlg, bg=biosim_theme('bg'))
        outer.pack(fill='both', expand=True, padx=10, pady=8)
        canv = tk.Canvas(outer, bg=biosim_theme('bg'), highlightthickness=0)
        sb = tk.Scrollbar(outer, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True)
        sb.pack(side='right', fill='y')
        inner = tk.Frame(canv, bg=biosim_theme('bg'))
        inner.bind('<Configure>',
                    lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=inner, anchor='nw', width=860)

        # Temporarily swap inner for our renderer
        old_inner = self.align_results_inner
        self.align_results_inner = inner
        try:
            self._render_alignment_result(hit['alignment'], local=True,
                                              name1='Query', name2=hit['name'][:6])
        finally:
            self.align_results_inner = old_inner

        tk.Button(dlg, text='Close',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 10), cursor='hand2', pady=4,
                   command=dlg.destroy
                   ).pack(fill='x', padx=14, pady=8)

    def _align_show_error(self, msg):
        self._align_clear_results()
        tk.Label(self.align_results_inner, text='' + msg,
                  bg=biosim_theme('bg'), fg='#cf222e',
                  font=('Arial', 11, 'bold'),
                  padx=20, pady=20).pack(anchor='center')

    def _blast_score_color(self, identity):
        if identity >= 80: return '#1a7f37'
        if identity >= 50: return '#9a6700'
        if identity >= 30: return '#bc4c00'
        return '#cf222e'

    def _build_tab_databases(self):
        """Mol-4: UniProt + NCBI Entrez + protein motif detection."""
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=self.T['tab_database'])
        self.tab_databases = tab

        # ─── Top: fetch panel ───
        ctrl = tk.LabelFrame(tab,
                                text='  ' + ('Загрузка из БД' if self.lang == 'ru'
                                                  else 'Database Fetch') + '  ',
                                bg=biosim_theme('bg'), fg='#8250df',
                                font=('Arial', 11, 'bold'), relief='flat')
        ctrl.pack(fill='x', padx=4, pady=4)

        # Mode buttons
        mode_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        mode_frame.pack(fill='x', padx=8, pady=4)
        tk.Label(mode_frame, text=('База:' if self.lang == 'ru' else 'Database:'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 10, 'bold')).pack(side='left', padx=(0, 8))

        self.db_source_var = tk.StringVar(value='uniprot')
        for val, lbl, col in [
            ('uniprot',     'UniProt (protein)',                      '#0969da'),
            ('ncbi_nuc',    'NCBI nucleotide',                         '#1a7f37'),
            ('ncbi_prot',   'NCBI protein',                            '#8250df'),
            ('motifs',      '' + ('Поиск мотивов' if self.lang == 'ru' else 'Find motifs'), '#9a6700'),
        ]:
            tk.Radiobutton(mode_frame, text=lbl, value=val,
                            variable=self.db_source_var,
                            bg=biosim_theme('bg'), fg=col,
                            selectcolor='#eef1f5',
                            activebackground='#f5f7fa',
                            font=('Arial', 9, 'bold'),
                            command=self._switch_db_mode
                            ).pack(side='left', padx=4)

        # Input row
        input_row = tk.Frame(ctrl, bg=biosim_theme('bg'))
        input_row.pack(fill='x', padx=8, pady=6)

        tk.Label(input_row, text=('ID / Accession:' if self.lang == 'ru'
                                       else 'ID / Accession:'),
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 10, 'bold')).pack(side='left')

        self.db_query_var = tk.StringVar()
        tk.Entry(input_row, textvariable=self.db_query_var,
                  bg=biosim_theme('panel2'), fg=biosim_theme('text'), insertbackground='#1c2128',
                  relief='flat', font=('Courier', 11),
                  highlightthickness=1, highlightbackground=biosim_theme('border'),
                  width=24
                  ).pack(side='left', padx=8, ipady=3)

        tk.Button(input_row, text='Fetch',
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4,
                   command=self._db_fetch
                   ).pack(side='left', padx=4)

        tk.Button(input_row, text='Search',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4,
                   command=self._db_search
                   ).pack(side='left', padx=4)

        tk.Button(input_row, text='Current sequence',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=10, pady=4,
                   command=self._db_use_current
                   ).pack(side='right', padx=4)

        # Hint
        hint_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        hint_frame.pack(fill='x', padx=8, pady=2)
        self.db_hint_lbl = tk.Label(hint_frame,
                                         text='',
                                         bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                                         font=('Arial', 9, 'italic'),
                                         wraplength=900, justify='left')
        self.db_hint_lbl.pack(anchor='w')

        # ─── Results area ───
        results_frame = tk.LabelFrame(tab,
                                          text='  ' + ('Результаты' if self.lang == 'ru'
                                                           else 'Results') + '  ',
                                          bg=biosim_theme('bg'), fg='#1a7f37',
                                          font=('Arial', 11, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=4, pady=4)

        canv = tk.Canvas(results_frame, bg=biosim_theme('bg'), highlightthickness=0)
        sb = tk.Scrollbar(results_frame, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        sb.pack(side='right', fill='y')
        self.db_results_inner = tk.Frame(canv, bg=biosim_theme('bg'))
        self.db_results_inner.bind('<Configure>',
                                       lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=self.db_results_inner,
                             anchor='nw', width=1100)

        self._switch_db_mode()  # set initial hint

    def _switch_db_mode(self):
        mode = self.db_source_var.get()
        hints = {
            'uniprot':   'Введи UniProt ID (например P04637 = TP53)' if self.lang == 'ru'
                          else 'Enter UniProt ID (e.g. P04637 = TP53)',
            'ncbi_nuc':  'Введи accession (например NM_001126112 или 12345) для нуклеотидов' if self.lang == 'ru'
                          else 'Enter accession (e.g. NM_001126112 or 12345) for nucleotides',
            'ncbi_prot': 'Введи accession белка (например NP_000537)' if self.lang == 'ru'
                          else 'Enter protein accession (e.g. NP_000537)',
            'motifs':    'Жми "Current sequence" чтобы искать мотивы в текущей последовательности' if self.lang == 'ru'
                          else 'Click "Current sequence" to find motifs in current sequence',
        }
        self.db_hint_lbl.config(text=hints.get(mode, ''))
        self._db_clear_results()

    def _db_clear_results(self):
        for w in self.db_results_inner.winfo_children():
            w.destroy()

    def _db_use_current(self):
        mode = self.db_source_var.get()
        if mode == 'motifs':
            # Find motifs in current sequence
            self._do_motif_search()
        else:
            # Load current sequence ID-style (placeholder, just show info)
            self.db_query_var.set(self.sequence_name or '')

    def _db_fetch(self):
        """Fetch sequence by accession from selected DB."""
        mode = self.db_source_var.get()
        acc = self.db_query_var.get().strip()
        if not acc:
            self._db_show_error('Enter an accession ID')
            return

        if mode == 'motifs':
            # If user pressed Fetch in motif mode, search motifs
            self._do_motif_search()
            return

        self.status_lbl.config(text=('Загружаю...' if self.lang == 'ru'
                                          else 'Fetching...'),
                                    fg='#0969da')
        self.update()
        self._db_clear_results()

        # Show loading
        loading = tk.Label(self.db_results_inner,
                              text='⏳ ' + ('Запрос к серверу...' if self.lang == 'ru'
                                              else 'Querying server...'),
                              bg=biosim_theme('bg'), fg='#9a6700',
                              font=('Arial', 11, 'italic'), pady=10)
        loading.pack()
        self.update()

        try:
            if mode == 'uniprot':
                result = fetch_uniprot_protein(acc)
                seq_type = 'protein'
            elif mode == 'ncbi_nuc':
                result = fetch_ncbi_sequence(acc, db='nucleotide')
                seq_type = 'dna'
            elif mode == 'ncbi_prot':
                result = fetch_ncbi_sequence(acc, db='protein')
                seq_type = 'protein'
            else:
                result = None
                seq_type = 'dna'
        except Exception as e:
            self._db_show_error('Fetch error: ' + str(e)[:80])
            biosim_log('error', 'db_fetch: ' + str(e))
            return

        loading.destroy()

        if result is None:
            self._db_show_error(
                ('Не удалось загрузить. Проверь ID и подключение к интернету.'
                  if self.lang == 'ru' else
                  'Could not fetch. Check ID and internet connection.'))
            return

        # Show result
        self._render_db_result(result, seq_type)
        self.status_lbl.config(
            text='Loaded ' + result['name'] + ' · ' + str(len(result['sequence'])) + ' bp/aa',
            fg='#1a7f37')

    def _db_search(self):
        """Search NCBI by query."""
        mode = self.db_source_var.get()
        query = self.db_query_var.get().strip()
        if not query:
            self._db_show_error('Enter a search query')
            return

        if mode not in ('ncbi_nuc', 'ncbi_prot'):
            self._db_show_error('Search is only for NCBI databases')
            return

        self.status_lbl.config(text=('Поиск...' if self.lang == 'ru'
                                          else 'Searching...'),
                                    fg='#0969da')
        self.update()
        self._db_clear_results()

        loading = tk.Label(self.db_results_inner,
                              text='⏳ ' + ('Поиск в NCBI...' if self.lang == 'ru'
                                              else 'Searching NCBI...'),
                              bg=biosim_theme('bg'), fg='#9a6700',
                              font=('Arial', 11, 'italic'), pady=10)
        loading.pack()
        self.update()

        try:
            results = search_ncbi_gene(query, max_results=15)
        except Exception as e:
            biosim_log('warn', '_db_search: {}'.format(e))
            self._db_show_error('Search error: ' + str(e)[:80])
            return

        loading.destroy()

        if not results:
            self._db_show_error('No results found')
            return

        # Header
        tk.Label(self.db_results_inner,
                  text='' + 'Найдено: ' + str(len(results)) + ' (NCBI)',
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 13, 'bold'), pady=4).pack(anchor='w', padx=12)

        # Each result as a row
        for r in results:
            card = tk.Frame(self.db_results_inner, bg=biosim_theme('panel'),
                              highlightthickness=1,
                              highlightbackground=biosim_theme('border'))
            card.pack(fill='x', padx=10, pady=2)
            head = tk.Frame(card, bg=biosim_theme('panel'))
            head.pack(fill='x', padx=10, pady=4)
            tk.Label(head, text=r['id'],
                      bg=biosim_theme('panel'), fg='#1a7f37',
                      font=('Courier', 11, 'bold'), width=14, anchor='w'
                      ).pack(side='left')
            tk.Label(head, text=r.get('title', '')[:80],
                      bg=biosim_theme('panel'), fg=biosim_theme('text'),
                      font=('Arial', 10), anchor='w'
                      ).pack(side='left', padx=8)

            def make_fetch(acc, db=mode):
                def fetch_it():
                    self.db_query_var.set(acc)
                    self._db_fetch()
                return fetch_it

            tk.Button(head, text='Fetch',
                       bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                       font=('Arial', 8, 'bold'), cursor='hand2',
                       padx=8, pady=1,
                       command=make_fetch(r['id'])
                       ).pack(side='right', padx=4)

    def _render_db_result(self, result, seq_type):
        """Display a fetched sequence with metadata."""
        # Info card
        card = tk.LabelFrame(self.db_results_inner,
                                text='  ' + result.get('name', '?') + '  ',
                                bg=biosim_theme('bg'), fg='#1a7f37',
                                font=('Arial', 11, 'bold'), relief='flat')
        card.pack(fill='x', padx=10, pady=4)

        # Metadata
        meta = tk.Frame(card, bg=biosim_theme('bg'))
        meta.pack(fill='x', padx=10, pady=6)

        info_items = [
            ('Source:', result.get('source', '?')),
            ('Accession:', result.get('accession', '?')),
            ('Length:', str(len(result['sequence'])) + (' aa' if seq_type == 'protein' else ' bp')),
        ]
        if result.get('organism'):
            info_items.append(('Organism:', result['organism']))
        if result.get('description'):
            info_items.append(('Description:', result['description'][:120]))

        for i, (k, v) in enumerate(info_items):
            tk.Label(meta, text=k, bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Arial', 9, 'bold'), anchor='e', width=12
                      ).grid(row=i, column=0, sticky='e', padx=4, pady=1)
            tk.Label(meta, text=v, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                      font=('Courier', 9), anchor='w', wraplength=850, justify='left'
                      ).grid(row=i, column=1, sticky='w', padx=4, pady=1)

        # Sequence preview
        preview = tk.LabelFrame(self.db_results_inner,
                                    text='  Sequence preview  ',
                                    bg=biosim_theme('bg'), fg='#0969da',
                                    font=('Arial', 10, 'bold'), relief='flat')
        preview.pack(fill='x', padx=10, pady=4)

        prev_text = tk.Text(preview, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                font=('Courier', 9), relief='flat',
                                wrap='word', height=6, padx=10, pady=8)
        prev_text.pack(fill='x', padx=4, pady=4)
        # Show first 480 chars
        seq = result['sequence']
        prev = format_sequence_blocks(seq[:480], block_size=10, line_size=60)
        prev_text.insert('1.0', prev)
        if len(seq) > 480:
            prev_text.insert('end', '\n... (' + str(len(seq) - 480) + ' more)')
        prev_text.config(state='disabled')

        # Action buttons
        btn_row = tk.Frame(self.db_results_inner, bg=biosim_theme('bg'))
        btn_row.pack(fill='x', padx=10, pady=4)

        def load_into_editor():
            self.current_sequence = seq
            self.sequence_name = result.get('name', 'fetched')
            self.sequence_type = seq_type
            if hasattr(self, 'seq_text'):
                self.seq_text.delete('1.0', 'end')
                self.seq_text.insert('1.0', seq)
                if hasattr(self, 'seq_name_var'):
                    self.seq_name_var.set(result.get('name', 'fetched'))
                if hasattr(self, 'seq_type_var'):
                    self.seq_type_var.set(seq_type)
                self._refresh_seq_analysis()
                self._refresh_seq_highlights()
                self.notebook.select(self.tab_sequence)
                self.status_lbl.config(
                    text='Loaded ' + str(len(seq)) + ' into Sequence Editor',
                    fg='#1a7f37')

        def find_motifs_here():
            if seq_type == 'protein':
                self._do_motif_search(seq=seq, name=result.get('name', '?'))
            else:
                self.status_lbl.config(text='Motif search needs protein',
                                            fg='#9a6700')

        def copy_seq():
            self.clipboard_clear()
            self.clipboard_append(seq)
            self.status_lbl.config(text='Sequence copied', fg='#1a7f37')

        tk.Button(btn_row, text='Load into editor',
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=4, command=load_into_editor
                   ).pack(side='left', padx=2)
        if seq_type == 'protein':
            tk.Button(btn_row, text='Find motifs',
                       bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                       font=('Arial', 10, 'bold'), cursor='hand2',
                       padx=12, pady=4, command=find_motifs_here
                       ).pack(side='left', padx=2)
        tk.Button(btn_row, text='Copy',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=12, pady=4, command=copy_seq
                   ).pack(side='left', padx=2)

    def _do_motif_search(self, seq=None, name=None):
        """Find protein motifs/domains in current or given sequence."""
        if seq is None:
            seq = self.current_sequence
            name = self.sequence_name or 'current'
        if not seq:
            self._db_show_error('No sequence loaded')
            return
        if self.sequence_type != 'protein' and not all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in seq[:50].upper()):
            # Auto-detect if it looks like protein
            non_dna = sum(1 for c in seq[:100].upper()
                            if c not in 'ATGCNU')
            if non_dna < 10:
                self._db_show_error('This appears to be DNA. Translate first or use protein sequence.')
                return

        motifs = find_protein_motifs(seq)
        self._db_clear_results()

        # Header
        tk.Label(self.db_results_inner,
                  text='' + 'Найдено мотивов: ' + str(len(motifs)) + ' в ' + str(name),
                  bg=biosim_theme('bg'), fg='#9a6700',
                  font=('Arial', 13, 'bold'), pady=4).pack(anchor='w', padx=12)

        if not motifs:
            tk.Label(self.db_results_inner,
                      text=('Известных мотивов не найдено. Попробуй на больших белках (TP53, GFP).'
                            if self.lang == 'ru' else
                            'No known motifs found. Try on bigger proteins (TP53, GFP).'),
                      bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Arial', 10, 'italic'), pady=10).pack()
            return

        # Group by motif
        by_motif = {}
        for m in motifs:
            by_motif.setdefault(m['motif'], []).append(m)

        for motif_name, occurrences in by_motif.items():
            card = tk.Frame(self.db_results_inner, bg=biosim_theme('panel'),
                              highlightthickness=2,
                              highlightbackground='#9a6700')
            card.pack(fill='x', padx=10, pady=4)

            head = tk.Frame(card, bg=biosim_theme('panel'))
            head.pack(fill='x', padx=10, pady=(8, 2))
            tk.Label(head, text='' + motif_name,
                      bg=biosim_theme('panel'), fg='#9a6700',
                      font=('Arial', 12, 'bold')).pack(side='left')
            tk.Label(head, text='  × {}'.format(len(occurrences)),
                      bg=biosim_theme('panel'), fg='#1a7f37',
                      font=('Courier', 11, 'bold')).pack(side='left', padx=8)

            tk.Label(card, text=occurrences[0]['description'],
                      bg=biosim_theme('panel'), fg=biosim_theme('text'),
                      font=('Arial', 9, 'italic'),
                      wraplength=900, justify='left',
                      ).pack(anchor='w', padx=14, pady=2)

            # Positions
            for m in occurrences[:5]:
                pos_text = '  position {}: '.format(m['position'] + 1) + m['match']
                tk.Label(card, text=pos_text,
                          bg=biosim_theme('panel'), fg='#0550ae',
                          font=('Courier', 10)).pack(anchor='w', padx=14, pady=0)
            if len(occurrences) > 5:
                tk.Label(card, text='  ... and ' + str(len(occurrences) - 5) + ' more',
                          bg=biosim_theme('panel'), fg=biosim_theme('text_dim'),
                          font=('Courier', 9, 'italic')
                          ).pack(anchor='w', padx=14)

            # Padding
            tk.Frame(card, bg=biosim_theme('panel'), height=6).pack()

    def _db_show_error(self, msg):
        self._db_clear_results()
        tk.Label(self.db_results_inner, text='' + msg,
                  bg=biosim_theme('bg'), fg='#cf222e',
                  font=('Arial', 11, 'bold'),
                  padx=20, pady=20).pack(anchor='center')

    def _build_tab_structure(self):
        """Mol-5: RNA structure, codon optimization, protein analysis."""
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=self.T['tab_structure'])
        self.tab_structure = tab

        # ─── Top controls ───
        ctrl = tk.LabelFrame(tab,
                                text='  ' + ('Структура и оптимизация' if self.lang == 'ru'
                                                  else 'Structure & Optimization') + '  ',
                                bg=biosim_theme('bg'), fg='#0969da',
                                font=('Arial', 11, 'bold'), relief='flat')
        ctrl.pack(fill='x', padx=4, pady=4)

        # Mode buttons
        mode_frame = tk.Frame(ctrl, bg=biosim_theme('bg'))
        mode_frame.pack(fill='x', padx=8, pady=4)
        tk.Label(mode_frame, text=('Анализ:' if self.lang == 'ru' else 'Analysis:'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 10, 'bold')).pack(side='left', padx=(0, 8))

        self.struct_mode_var = tk.StringVar(value='rna')
        for val, lbl, col in [
            ('rna',     '' + ('РНК структура (Nussinov)' if self.lang == 'ru' else 'RNA structure (Nussinov)'), '#1a7f37'),
            ('codon',   '' + ('Кодон оптимизация' if self.lang == 'ru' else 'Codon optimization'), '#9a6700'),
            ('protein', '' + ('Анализ белка' if self.lang == 'ru' else 'Protein analysis'), '#0969da'),
            ('esmfold', 'ESMFold (3D)', '#8250df'),
        ]:
            tk.Radiobutton(mode_frame, text=lbl, value=val,
                            variable=self.struct_mode_var,
                            bg=biosim_theme('bg'), fg=col,
                            selectcolor='#eef1f5',
                            activebackground='#f5f7fa',
                            font=('Arial', 9, 'bold'),
                            command=self._switch_struct_mode
                            ).pack(side='left', padx=4)

        # Action row
        action_row = tk.Frame(ctrl, bg=biosim_theme('bg'))
        action_row.pack(fill='x', padx=8, pady=6)

        # Codon optimization: organism selector
        self.struct_organism_var = tk.StringVar(value='ecoli')
        self.struct_organism_frame = tk.Frame(action_row, bg=biosim_theme('bg'))
        tk.Label(self.struct_organism_frame,
                  text=('Организм:' if self.lang == 'ru' else 'Organism:'),
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 10, 'bold')).pack(side='left', padx=(0, 4))
        org_dd = ttk.Combobox(self.struct_organism_frame,
                                 textvariable=self.struct_organism_var,
                                 values=['ecoli', 'yeast', 'human'],
                                 state='readonly', width=10,
                                 font=('Arial', 10))
        org_dd.pack(side='left', padx=2)
        # Initially packed if mode is codon
        self.struct_organism_frame.pack(side='left', padx=4)

        tk.Button(action_row,
                   text='▶ ' + ('Запустить' if self.lang == 'ru' else 'Run'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=16, pady=4,
                   command=self._run_structure
                   ).pack(side='right', padx=4)

        tk.Button(action_row,
                   text=('Использовать текущую' if self.lang == 'ru'
                                  else 'Use current sequence'),
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=10, pady=4,
                   command=self._struct_use_current
                   ).pack(side='right', padx=4)

        # Hint label
        self.struct_hint_lbl = tk.Label(ctrl, text='',
                                            bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                                            font=('Arial', 9, 'italic'),
                                            wraplength=900, justify='left')
        self.struct_hint_lbl.pack(anchor='w', padx=10, pady=2)

        # ─── Results area ───
        results_frame = tk.LabelFrame(tab,
                                          text='  ' + ('Результаты' if self.lang == 'ru'
                                                           else 'Results') + '  ',
                                          bg=biosim_theme('bg'), fg='#1a7f37',
                                          font=('Arial', 11, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=4, pady=4)

        canv = tk.Canvas(results_frame, bg=biosim_theme('bg'), highlightthickness=0)
        sb = tk.Scrollbar(results_frame, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        sb.pack(side='right', fill='y')
        self.struct_results_inner = tk.Frame(canv, bg=biosim_theme('bg'))
        self.struct_results_inner.bind('<Configure>',
                                            lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=self.struct_results_inner,
                             anchor='nw', width=1100)

        self._switch_struct_mode()

    def _switch_struct_mode(self):
        mode = self.struct_mode_var.get()
        hints = {
            'rna':     ('Предсказание вторичной структуры РНК (Nussinov, max pairs). Лимит 500 nt.'
                         if self.lang == 'ru' else
                         'RNA secondary structure prediction (Nussinov, max base pairs). Limit 500 nt.'),
            'codon':   ('Оптимизация кодонов под организм для усиления экспрессии.'
                         if self.lang == 'ru' else
                         'Codon optimization for chosen organism (boost expression).'),
            'protein': ('Анализ белка: pI, MW, extinction coefficient, гидрофобность.'
                         if self.lang == 'ru' else
                         'Protein analysis: pI, MW, extinction coefficient, hydrophobicity profile.'),
            'esmfold': ('Предсказание 3D структуры белка через ESMFold (требует PyTorch + ESM модель).'
                         if self.lang == 'ru' else
                         'Predict 3D protein structure via ESMFold (requires PyTorch + ESM model).'),
        }
        self.struct_hint_lbl.config(text=hints.get(mode, ''))
        # Show organism only for codon mode
        if mode == 'codon':
            self.struct_organism_frame.pack(side='left', padx=4)
        else:
            self.struct_organism_frame.pack_forget()
        self._struct_clear_results()

    def _struct_clear_results(self):
        for w in self.struct_results_inner.winfo_children():
            w.destroy()

    def _struct_use_current(self):
        if not self.current_sequence:
            self._struct_show_error('No sequence loaded')
            return
        self._run_structure()

    def _run_structure(self):
        mode = self.struct_mode_var.get()
        if not self.current_sequence:
            self._struct_show_error('Load a sequence first (in Sequence Editor tab)'
                                        if self.lang == 'en' else
                                        'Сначала загрузи последовательность (в редакторе)')
            return

        if mode == 'rna':
            self._do_rna_structure()
        elif mode == 'codon':
            self._do_codon_optimize()
        elif mode == 'protein':
            self._do_protein_analysis()
        elif mode == 'esmfold':
            self._do_esmfold()

    def _do_rna_structure(self):
        """Predict and display RNA secondary structure."""
        seq = self.current_sequence
        if self.sequence_type == 'protein':
            self._struct_show_error('Protein sequence — need DNA/RNA')
            return

        self.status_lbl.config(text='Predicting RNA structure...' if self.lang == 'en'
                                    else 'Предсказание структуры РНК...',
                                    fg='#0969da')
        self.update()

        try:
            result = predict_rna_structure_nussinov(seq, min_loop=3)
        except Exception as e:
            biosim_log('warn', '_do_rna_structure: {}'.format(e))
            self._struct_show_error('Error: ' + str(e)[:80])
            return

        if result is None:
            self._struct_show_error('Sequence too long (max 500 nt for Nussinov)')
            return

        self._struct_clear_results()

        # Header card
        head = tk.LabelFrame(self.struct_results_inner,
                                text='  ' + ('РНК вторичная структура' if self.lang == 'ru'
                                                  else 'RNA secondary structure') + '  ',
                                bg=biosim_theme('bg'), fg='#1a7f37',
                                font=('Arial', 11, 'bold'), relief='flat')
        head.pack(fill='x', padx=10, pady=4)

        stats = tk.Frame(head, bg=biosim_theme('bg'))
        stats.pack(fill='x', padx=10, pady=6)

        items = [
            ('Length:', '{} nt'.format(len(result['sequence']))),
            ('Pairs:', '{}'.format(result['score'])),
            ('Pair density:',
              '{:.2f}%'.format(100 * 2 * result['score'] / max(1, len(result['sequence'])))),
            ('ΔG estimate:', '{:.1f} kcal/mol'.format(result['dG_estimate'])),
        ]
        for i, (k, v) in enumerate(items):
            tk.Label(stats, text=k, bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Arial', 10), anchor='e', width=14
                      ).grid(row=i, column=0, sticky='e', padx=4, pady=1)
            tk.Label(stats, text=v, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                      font=('Courier', 10, 'bold'), anchor='w'
                      ).grid(row=i, column=1, sticky='w', padx=4, pady=1)

        # Dot-bracket display
        db_frame = tk.LabelFrame(self.struct_results_inner,
                                    text='  Dot-bracket notation  ',
                                    bg=biosim_theme('bg'), fg='#0550ae',
                                    font=('Arial', 10, 'bold'), relief='flat')
        db_frame.pack(fill='x', padx=10, pady=4)

        db_text_frame = tk.Frame(db_frame, bg=biosim_theme('panel'))
        db_text_frame.pack(fill='x', padx=4, pady=4)
        db_sb = tk.Scrollbar(db_text_frame, orient='horizontal')
        db_sb.pack(side='bottom', fill='x')
        db_text = tk.Text(db_text_frame, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                              font=('Courier', 11), relief='flat',
                              wrap='none', height=6, padx=10, pady=8,
                              xscrollcommand=db_sb.set)
        db_text.pack(side='top', fill='x', expand=True)
        db_sb.config(command=db_text.xview)

        db_text.tag_config('paired', foreground='#1a7f37',
                              font=('Courier', 11, 'bold'))
        db_text.tag_config('unpaired', foreground=biosim_theme('text_dim'),
                              font=('Courier', 11))
        db_text.tag_config('nucl_paired', foreground='#9a6700',
                              font=('Courier', 11, 'bold'))
        db_text.tag_config('nucl_unpaired', foreground=biosim_theme('text'),
                              font=('Courier', 11))

        chunk = 60
        for s in range(0, len(result['sequence']), chunk):
            seq_chunk = result['sequence'][s:s+chunk]
            db_chunk = result['dot_bracket'][s:s+chunk]
            # Position
            db_text.insert('end', '{:>5} '.format(s + 1), 'unpaired')
            # Sequence
            for i, ch in enumerate(seq_chunk):
                tag = 'nucl_paired' if db_chunk[i] != '.' else 'nucl_unpaired'
                db_text.insert('end', ch, tag)
            db_text.insert('end', '\n')
            db_text.insert('end', '      ')
            for ch in db_chunk:
                tag = 'paired' if ch != '.' else 'unpaired'
                db_text.insert('end', ch, tag)
            db_text.insert('end', '\n\n')
        db_text.config(state='disabled')

        # Visual structure (arc diagram)
        if result['pairs']:
            arc_frame = tk.LabelFrame(self.struct_results_inner,
                                          text='  ' + ('Дуговая диаграмма' if self.lang == 'ru'
                                                            else 'Arc diagram') + '  ',
                                          bg=biosim_theme('bg'), fg='#8250df',
                                          font=('Arial', 10, 'bold'), relief='flat')
            arc_frame.pack(fill='x', padx=10, pady=4)

            canv_w = 1000
            canv_h = 250
            c = tk.Canvas(arc_frame, bg=biosim_theme('panel'), width=canv_w, height=canv_h,
                            highlightthickness=0)
            c.pack(padx=4, pady=4)
            
            seq_str = result['sequence']
            n = len(seq_str)
            margin = 40
            plot_w = canv_w - 2 * margin
            baseline = canv_h - 40
            
            # Draw sequence as baseline
            for i in range(n):
                x = margin + plot_w * i / max(1, n - 1)
                # Nucleotide color
                nt = seq_str[i]
                col_map = {'A': '#1a7f37', 'U': '#cf222e',
                            'G': '#9a6700', 'C': '#0969da'}
                color = col_map.get(nt, '#57606a')
                c.create_oval(x - 2, baseline - 2, x + 2, baseline + 2,
                                fill=color, outline='')
                # Labels every 10
                if i % 10 == 0:
                    c.create_text(x, baseline + 14,
                                    text=str(i + 1), fill='#57606a',
                                    font=('Courier', 7))
            
            # Draw arcs for pairs
            import math as _m
            for (i, j) in result['pairs']:
                x1 = margin + plot_w * i / max(1, n - 1)
                x2 = margin + plot_w * j / max(1, n - 1)
                # Arc above baseline
                arc_height = (x2 - x1) * 0.5
                cx = (x1 + x2) / 2
                cy = baseline
                # Bounding box for arc
                r = (x2 - x1) / 2
                c.create_arc(cx - r, cy - r, cx + r, cy + r,
                                start=0, extent=180,
                                style='arc', outline='#8250df', width=1)

        # Copy buttons
        btn_row = tk.Frame(self.struct_results_inner, bg=biosim_theme('bg'))
        btn_row.pack(fill='x', padx=10, pady=4)

        def copy_dotbracket():
            text = '> ' + (self.sequence_name or 'sequence') + '\n'
            text += result['sequence'] + '\n'
            text += result['dot_bracket'] + '\n'
            self.clipboard_clear()
            self.clipboard_append(text)
            self.status_lbl.config(text='Dot-bracket copied', fg='#1a7f37')

        tk.Button(btn_row, text='Copy (Vienna format)',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=10, pady=3, command=copy_dotbracket
                   ).pack(side='left', padx=2)

    def _do_codon_optimize(self):
        """Codon optimization for selected organism."""
        seq = self.current_sequence
        organism = self.struct_organism_var.get()

        if self.sequence_type == 'protein':
            # Direct optimization
            optimized = codon_optimize(seq, organism)
            original_cai = None
        else:
            # DNA — first compute CAI, then translate + optimize
            original_cai = cai(seq, organism)
            protein = translate_dna(seq, frame=0).rstrip('*')
            if not protein:
                self._struct_show_error('Cannot translate sequence')
                return
            optimized = codon_optimize(protein, organism)

        new_cai = cai(optimized, organism)

        self._struct_clear_results()

        # Stats
        head = tk.LabelFrame(self.struct_results_inner,
                                text='  ' + ('Оптимизация кодонов' if self.lang == 'ru'
                                                  else 'Codon Optimization') + '  ',
                                bg=biosim_theme('bg'), fg='#9a6700',
                                font=('Arial', 11, 'bold'), relief='flat')
        head.pack(fill='x', padx=10, pady=4)

        stats = tk.Frame(head, bg=biosim_theme('bg'))
        stats.pack(fill='x', padx=10, pady=6)

        items = [
            ('Organism:', organism),
            ('Optimized length:', '{} nt = {} aa'.format(
                len(optimized), len(optimized) // 3)),
            ('CAI (new):', '{:.3f}'.format(new_cai)),
        ]
        if original_cai is not None:
            items.insert(2, ('CAI (original):', '{:.3f}'.format(original_cai)))
            items.append(('CAI improvement:',
                            '+{:.3f}'.format(new_cai - original_cai)))

        for i, (k, v) in enumerate(items):
            tk.Label(stats, text=k, bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Arial', 10), anchor='e', width=18
                      ).grid(row=i, column=0, sticky='e', padx=4, pady=1)
            tk.Label(stats, text=v, bg=biosim_theme('bg'), fg='#1a7f37',
                      font=('Courier', 10, 'bold'), anchor='w'
                      ).grid(row=i, column=1, sticky='w', padx=4, pady=1)

        # Optimized sequence
        opt_frame = tk.LabelFrame(self.struct_results_inner,
                                      text='  Optimized DNA  ',
                                      bg=biosim_theme('bg'), fg='#1a7f37',
                                      font=('Arial', 10, 'bold'), relief='flat')
        opt_frame.pack(fill='x', padx=10, pady=4)

        opt_text = tk.Text(opt_frame, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                              font=('Courier', 10), relief='flat',
                              wrap='word', height=8, padx=10, pady=8)
        opt_text.pack(fill='x', padx=4, pady=4)
        opt_text.insert('1.0', format_sequence_blocks(optimized,
                                                          block_size=10, line_size=60))
        opt_text.config(state='disabled')

        # Buttons
        btn_row = tk.Frame(self.struct_results_inner, bg=biosim_theme('bg'))
        btn_row.pack(fill='x', padx=10, pady=4)

        def copy_opt():
            self.clipboard_clear()
            self.clipboard_append(optimized)
            self.status_lbl.config(text='Optimized DNA copied', fg='#1a7f37')

        def load_opt():
            self.current_sequence = optimized
            self.sequence_type = 'dna'
            self.sequence_name = (self.sequence_name or 'gene') + '_opt_' + organism
            if hasattr(self, 'seq_text'):
                self.seq_text.delete('1.0', 'end')
                self.seq_text.insert('1.0', optimized)
                if hasattr(self, 'seq_name_var'):
                    self.seq_name_var.set(self.sequence_name)
                if hasattr(self, 'seq_type_var'):
                    self.seq_type_var.set('dna')
                self._refresh_seq_analysis()
                self._refresh_seq_highlights()
                self.notebook.select(self.tab_sequence)
                self.status_lbl.config(text='Loaded optimized DNA',
                                            fg='#1a7f37')

        tk.Button(btn_row, text='Copy', bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'),
                   relief='solid', bd=1, font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=10, pady=3, command=copy_opt
                   ).pack(side='left', padx=2)
        tk.Button(btn_row, text='Load into editor',
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=10, pady=3, command=load_opt
                   ).pack(side='left', padx=2)

    def _do_protein_analysis(self):
        """Analyze protein properties."""
        seq = self.current_sequence
        if self.sequence_type != 'protein':
            # Try to translate
            if self.sequence_type in ('dna', 'rna'):
                seq = translate_dna(seq, frame=0).rstrip('*')
                if not seq:
                    self._struct_show_error('Cannot translate to protein')
                    return
            else:
                self._struct_show_error('Need protein sequence')
                return

        seq_clean = seq.replace('*', '')
        if not seq_clean:
            self._struct_show_error('No amino acids')
            return

        # Compute everything
        n = len(seq_clean)
        mw = molecular_weight_protein(seq_clean)
        pi = calculate_protein_pi(seq_clean)
        e_ox, e_red = calculate_extinction_coefficient(seq_clean)
        hp_profile = hydrophobicity_profile(seq_clean, window=9)

        # AA composition
        aa_counts = {}
        for aa in seq_clean:
            if aa in AA_PROPERTIES:
                aa_counts[aa] = aa_counts.get(aa, 0) + 1

        # Net charge at pH 7
        def net_charge_at_pH(pH):
            n_C = seq_clean.count('C')
            n_D = seq_clean.count('D')
            n_E = seq_clean.count('E')
            n_H = seq_clean.count('H')
            n_K = seq_clean.count('K')
            n_R = seq_clean.count('R')
            n_Y = seq_clean.count('Y')
            pos = (1 / (1 + 10 ** (pH - AA_PKA['N_TERM']))
                    + n_H * (1 / (1 + 10 ** (pH - AA_PKA['H'])))
                    + n_K * (1 / (1 + 10 ** (pH - AA_PKA['K'])))
                    + n_R * (1 / (1 + 10 ** (pH - AA_PKA['R']))))
            neg = (1 / (1 + 10 ** (AA_PKA['C_TERM'] - pH))
                    + n_D * (1 / (1 + 10 ** (AA_PKA['D'] - pH)))
                    + n_E * (1 / (1 + 10 ** (AA_PKA['E'] - pH)))
                    + n_C * (1 / (1 + 10 ** (AA_PKA['C'] - pH)))
                    + n_Y * (1 / (1 + 10 ** (AA_PKA['Y'] - pH))))
            return pos - neg

        charge_7 = net_charge_at_pH(7.0)
        charge_5 = net_charge_at_pH(5.0)

        self._struct_clear_results()

        # Main info card
        head = tk.LabelFrame(self.struct_results_inner,
                                text='  ' + ('Анализ белка' if self.lang == 'ru'
                                                  else 'Protein Analysis') + '  ',
                                bg=biosim_theme('bg'), fg='#0969da',
                                font=('Arial', 11, 'bold'), relief='flat')
        head.pack(fill='x', padx=10, pady=4)

        stats = tk.Frame(head, bg=biosim_theme('bg'))
        stats.pack(fill='x', padx=10, pady=6)

        items = [
            ('Length:', '{} aa'.format(n)),
            ('Molecular weight:', '{:.1f} Da ({:.2f} kDa)'.format(mw, mw / 1000)),
            ('Isoelectric point (pI):', '{:.2f}'.format(pi)),
            ('Charge at pH 7:', '{:+.2f}'.format(charge_7)),
            ('Charge at pH 5:', '{:+.2f}'.format(charge_5)),
            ('ε280 (oxidized):', '{} M⁻¹·cm⁻¹'.format(e_ox)),
            ('ε280 (reduced):', '{} M⁻¹·cm⁻¹'.format(e_red)),
            ('A280 (1 mg/mL):', '{:.3f}'.format(e_red / mw) if mw > 0 else '?'),
        ]

        for i, (k, v) in enumerate(items):
            col = i // 4
            row = i % 4
            tk.Label(stats, text=k, bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Arial', 9), anchor='e', width=20
                      ).grid(row=row, column=col * 2, sticky='e', padx=4, pady=1)
            tk.Label(stats, text=v, bg=biosim_theme('bg'), fg='#1a7f37',
                      font=('Courier', 9, 'bold'), anchor='w', width=22
                      ).grid(row=row, column=col * 2 + 1, sticky='w', padx=4, pady=1)

        # AA composition table
        comp_frame = tk.LabelFrame(self.struct_results_inner,
                                       text='  ' + ('Состав аминокислот' if self.lang == 'ru'
                                                        else 'Amino acid composition') + '  ',
                                       bg=biosim_theme('bg'), fg='#1a7f37',
                                       font=('Arial', 10, 'bold'), relief='flat')
        comp_frame.pack(fill='x', padx=10, pady=4)

        # Bar chart
        comp_canvas = tk.Canvas(comp_frame, height=180, width=1000,
                                   bg=biosim_theme('panel'), highlightthickness=0)
        comp_canvas.pack(padx=4, pady=4)

        aa_list = sorted(aa_counts.keys(),
                          key=lambda a: -aa_counts[a])
        bar_w = 1000 / max(20, len(aa_list))
        max_count = max(aa_counts.values()) if aa_counts else 1
        for i, aa in enumerate(aa_list):
            count = aa_counts[aa]
            pct = 100 * count / n
            bar_h = int(130 * count / max_count)
            x = 20 + i * bar_w
            aa_type = AA_PROPERTIES.get(aa, {}).get('type', '')
            color = AA_COLORS.get(aa_type, '#57606a')
            comp_canvas.create_rectangle(
                x, 150 - bar_h, x + bar_w - 4, 150,
                fill=color, outline='')
            comp_canvas.create_text(x + bar_w / 2, 160,
                                       text=aa, fill='#1c2128',
                                       font=('Arial', 9, 'bold'))
            comp_canvas.create_text(x + bar_w / 2, 145 - bar_h,
                                       text='{:.1f}%'.format(pct),
                                       fill='#57606a', font=('Arial', 7),
                                       anchor='s')

        # Hydrophobicity profile
        if hp_profile:
            hp_frame = tk.LabelFrame(self.struct_results_inner,
                                         text='  Hydrophobicity (Kyte-Doolittle, window=9)  ',
                                         bg=biosim_theme('bg'), fg='#bc4c00',
                                         font=('Arial', 10, 'bold'), relief='flat')
            hp_frame.pack(fill='x', padx=10, pady=4)

            hp_canvas = tk.Canvas(hp_frame, height=150, width=1000,
                                      bg=biosim_theme('panel'), highlightthickness=0)
            hp_canvas.pack(padx=4, pady=4)

            # Find min/max for scaling
            min_h = min(p[1] for p in hp_profile)
            max_h = max(p[1] for p in hp_profile)
            range_h = max(max_h - min_h, 0.1)

            # Zero line
            zero_y = 75 + 60 * (0 - (min_h + max_h) / 2) / range_h
            hp_canvas.create_line(20, zero_y, 980, zero_y,
                                     fill='#d0d7de', dash=(2, 2))
            hp_canvas.create_text(15, zero_y, text='0', fill='#57606a',
                                     font=('Arial', 8), anchor='e')

            # Plot
            prev_x, prev_y = None, None
            for (pos, h) in hp_profile:
                x = 20 + 960 * pos / max(1, len(seq_clean) - 1)
                y = 75 - 60 * (h - (min_h + max_h) / 2) / range_h
                # Color by sign
                color = '#bc4c00' if h > 0 else '#0969da'
                if prev_x is not None:
                    hp_canvas.create_line(prev_x, prev_y, x, y,
                                              fill=color, width=2)
                prev_x, prev_y = x, y
            
            hp_canvas.create_text(15, 15, text=('Hydrophobic ' if self.lang == 'en'
                                                    else 'Гидрофобный '),
                                     fill='#bc4c00', font=('Arial', 8, 'italic'),
                                     anchor='w')
            hp_canvas.create_text(15, 135, text=('Hydrophilic ' if self.lang == 'en'
                                                     else 'Гидрофильный '),
                                     fill='#0969da', font=('Arial', 8, 'italic'),
                                     anchor='w')

    def _do_esmfold(self):
        """Try ESMFold prediction."""
        seq = self.current_sequence
        if self.sequence_type != 'protein':
            self._struct_show_error('Need protein sequence')
            return

        self._struct_clear_results()
        tk.Label(self.struct_results_inner,
                  text='⏳ ' + ('Проверка ESMFold...' if self.lang == 'ru'
                                  else 'Checking ESMFold availability...'),
                  bg=biosim_theme('bg'), fg='#9a6700',
                  font=('Arial', 11, 'italic'), pady=10).pack()
        self.update()

        result = try_esmfold(seq)

        self._struct_clear_results()

        card = tk.LabelFrame(self.struct_results_inner,
                                text='  ESMFold  ',
                                bg=biosim_theme('bg'), fg='#8250df',
                                font=('Arial', 11, 'bold'), relief='flat')
        card.pack(fill='x', padx=10, pady=4)

        if result is None or not result.get('success'):
            msg = (result or {}).get('message', 'ESMFold not available')
            tk.Label(card, text='' + msg,
                      bg=biosim_theme('bg'), fg='#9a6700',
                      font=('Arial', 11), wraplength=900,
                      justify='left', pady=10).pack(padx=14)

            # Suggest online alternatives
            tk.Label(card,
                      text=('Используй онлайн сервисы:' if self.lang == 'ru'
                            else 'Use online services:'),
                      bg=biosim_theme('bg'), fg='#1a7f37',
                      font=('Arial', 10, 'bold')
                      ).pack(anchor='w', padx=14, pady=(8, 4))
            
            services = [
                ('ESMAtlas', 'https://esmatlas.com/',
                  'Predict any sequence online (Meta AI)'),
                ('AlphaFold DB', 'https://alphafold.ebi.ac.uk/',
                  '200M+ pre-computed structures (DeepMind)'),
                ('ColabFold', 'https://colab.research.google.com/github/sokrypton/ColabFold',
                  'Free Google Colab notebooks'),
                ('SWISS-MODEL', 'https://swissmodel.expasy.org/',
                  'Homology modeling'),
            ]
            for name, url, desc in services:
                row = tk.Frame(card, bg=biosim_theme('bg'))
                row.pack(fill='x', padx=14, pady=2)
                tk.Label(row, text='• ' + name,
                          bg=biosim_theme('bg'), fg='#0969da',
                          font=('Arial', 10, 'bold'), anchor='w', width=15
                          ).pack(side='left')
                tk.Label(row, text=desc,
                          bg=biosim_theme('bg'), fg=biosim_theme('text'),
                          font=('Arial', 9), anchor='w'
                          ).pack(side='left', padx=4)

                def make_open(u):
                    def open_url():
                        try:
                            import webbrowser
                            webbrowser.open(u)
                            self.status_lbl.config(text='Opened ' + u,
                                                        fg='#1a7f37')
                        except Exception as _exc:
                            biosim_log('warn', 'open_url: {}'.format(_exc))
                    return open_url

                tk.Button(row, text='Open',
                           bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                           font=('Arial', 8, 'bold'), cursor='hand2',
                           padx=8, pady=1,
                           command=make_open(url)
                           ).pack(side='right', padx=2)

            # Sequence info
            tk.Label(card,
                      text=('Твоя последовательность готова для копирования:'
                            if self.lang == 'ru' else
                            'Your sequence ready to paste:'),
                      bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Arial', 9, 'italic')
                      ).pack(anchor='w', padx=14, pady=(10, 2))
            seq_text = tk.Text(card, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                  font=('Courier', 9), relief='flat',
                                  wrap='word', height=4, padx=10, pady=6)
            seq_text.pack(fill='x', padx=14, pady=4)
            seq_text.insert('1.0', seq)
            seq_text.config(state='disabled')

            tk.Button(card, text='Copy sequence',
                       bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                       font=('Arial', 10, 'bold'), cursor='hand2',
                       padx=12, pady=3,
                       command=lambda: (self.clipboard_clear(),
                                           self.clipboard_append(seq),
                                           self.status_lbl.config(
                                               text='Sequence copied',
                                               fg='#1a7f37'))
                       ).pack(anchor='w', padx=14, pady=(4, 10))
        else:
            # ESMFold actually succeeded (very unlikely in normal envs)
            tk.Label(card,
                      text='' + result.get('message', 'Predicted'),
                      bg=biosim_theme('bg'), fg='#1a7f37',
                      font=('Arial', 11, 'bold')).pack(padx=14, pady=10)
            if result.get('pdb_string'):
                pdb_text = tk.Text(card, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                      font=('Courier', 8), relief='flat',
                                      wrap='word', height=12, padx=10, pady=6)
                pdb_text.pack(fill='x', padx=14, pady=4)
                pdb_text.insert('1.0', result['pdb_string'][:5000])
                pdb_text.config(state='disabled')


    def _build_tab_crispr(self):
        """CRISPR/Cas9 guide RNA designer."""
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=self.T['tab_crispr'])
        self.tab_crispr = tab

        # ─── Top controls ───
        ctrl = tk.LabelFrame(tab,
                                text='  ' + ('Дизайн CRISPR/Cas9 гайдов'
                                                  if self.lang == 'ru' else
                                                  'CRISPR/Cas9 sgRNA Designer') + '  ',
                                bg=biosim_theme('bg'), fg='#1a7f37',
                                font=('Arial', 11, 'bold'), relief='flat')
        ctrl.pack(fill='x', padx=4, pady=4)

        params_grid = tk.Frame(ctrl, bg=biosim_theme('bg'))
        params_grid.pack(fill='x', padx=8, pady=6)

        tk.Label(params_grid, text=('Cas нуклеаза:' if self.lang == 'ru'
                                       else 'Cas variant:'),
                  bg=biosim_theme('bg'), fg='#0550ae',
                  font=('Arial', 10, 'bold')).grid(row=0, column=0, sticky='e', padx=4)
        self.crispr_cas_var = tk.StringVar(value='SpCas9')
        cas_dd = ttk.Combobox(params_grid, textvariable=self.crispr_cas_var,
                                 values=['SpCas9 (NGG, 20nt)',
                                          'SaCas9 (NNGRRT, 21nt)',
                                          'Cas12a/Cpf1 (TTTV, 23nt)'],
                                 state='readonly', width=24,
                                 font=('Arial', 10))
        cas_dd.pack_forget()
        cas_dd.grid(row=0, column=1, sticky='w', padx=4)

        tk.Label(params_grid, text=('Использовать текущую последовательность:'
                                       if self.lang == 'ru' else
                                       'Use current sequence:'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9, 'italic')).grid(row=1, column=0, columnspan=2,
                                                         sticky='w', padx=4, pady=(8, 2))

        # Sequence input
        seq_frame = tk.Frame(ctrl, bg=biosim_theme('panel'))
        seq_frame.pack(fill='x', padx=8, pady=4)
        sb_crispr = tk.Scrollbar(seq_frame)
        sb_crispr.pack(side='right', fill='y')
        self.crispr_seq_text = tk.Text(seq_frame, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                            insertbackground='#1c2128',
                                            relief='flat', font=('Courier', 9),
                                            height=4, padx=8, pady=4,
                                            yscrollcommand=sb_crispr.set)
        self.crispr_seq_text.pack(side='left', fill='x', expand=True)
        sb_crispr.config(command=self.crispr_seq_text.yview)

        btn_row = tk.Frame(ctrl, bg=biosim_theme('bg'))
        btn_row.pack(fill='x', padx=8, pady=4)

        def use_current():
            if self.current_sequence:
                self.crispr_seq_text.delete('1.0', 'end')
                self.crispr_seq_text.insert('1.0', self.current_sequence)

        tk.Button(btn_row, text=('Из редактора' if self.lang == 'ru'
                                            else 'From editor'),
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=10, pady=2, command=use_current
                   ).pack(side='left', padx=4)

        tk.Button(btn_row, text='▶ ' + ('Дизайн гайдов' if self.lang == 'ru'
                                            else 'Design guides'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=3, command=self._run_crispr
                   ).pack(side='right', padx=4)

        # ─── Results ───
        results_frame = tk.LabelFrame(tab,
                                          text='  ' + ('Дизайн гайдов' if self.lang == 'ru'
                                                            else 'Guide candidates') + '  ',
                                          bg=biosim_theme('bg'), fg='#1a7f37',
                                          font=('Arial', 11, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=4, pady=4)

        canv = tk.Canvas(results_frame, bg=biosim_theme('bg'), highlightthickness=0)
        sb = tk.Scrollbar(results_frame, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        sb.pack(side='right', fill='y')
        self.crispr_results_inner = tk.Frame(canv, bg=biosim_theme('bg'))
        self.crispr_results_inner.bind('<Configure>',
                                            lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=self.crispr_results_inner,
                             anchor='nw', width=1100)

        tk.Label(self.crispr_results_inner,
                  text=('Загрузи ДНК и жми "Дизайн гайдов" — программа найдёт sgRNA с PAM, GC%, оценкой.'
                        if self.lang == 'ru' else
                        'Load DNA and click "Design guides" — program finds sgRNAs with PAM, GC%, score.'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 11, 'italic'),
                  padx=20, pady=20).pack(anchor='center')

    def _run_crispr(self):
        """Design CRISPR guides for the loaded sequence."""
        seq = self.crispr_seq_text.get('1.0', 'end').strip().upper()
        seq = ''.join(c for c in seq if c in 'ATGCN')
        if len(seq) < 30:
            for w in self.crispr_results_inner.winfo_children():
                w.destroy()
            tk.Label(self.crispr_results_inner,
                      text=('Нужна последовательность >30 нт' if self.lang == 'ru'
                                      else 'Need sequence >30 nt'),
                      bg=biosim_theme('bg'), fg='#cf222e',
                      font=('Arial', 11, 'bold'),
                      padx=20, pady=20).pack(anchor='center')
            return

        cas = self.crispr_cas_var.get()
        if 'SaCas9' in cas:
            pam = 'NNGRRT'
            guide_len = 21
        elif 'Cas12a' in cas or 'Cpf1' in cas:
            pam = 'TTTV'
            guide_len = 23
        else:
            pam = 'NGG'
            guide_len = 20

        # For non-standard PAMs use NGG (simplification)
        guides = design_crispr_guides(seq, pam='NGG', length=guide_len,
                                          max_guides=20)

        # Clear and display
        for w in self.crispr_results_inner.winfo_children():
            w.destroy()

        tk.Label(self.crispr_results_inner,
                  text=('Найдено гайдов' if self.lang == 'ru' else 'Guides found') +
                        ': {}'.format(len(guides)),
                  bg=biosim_theme('bg'), fg='#1a7f37',
                  font=('Arial', 13, 'bold'), pady=4).pack(anchor='w', padx=12)

        if not guides:
            tk.Label(self.crispr_results_inner,
                      text=('Подходящих PAM сайтов не найдено' if self.lang == 'ru'
                            else 'No suitable PAM sites found'),
                      bg=biosim_theme('bg'), fg='#9a6700',
                      font=('Arial', 11, 'italic'),
                      padx=20, pady=20).pack(anchor='center')
            return

        # Header table
        header = tk.Frame(self.crispr_results_inner, bg=biosim_theme('panel'))
        header.pack(fill='x', padx=10, pady=(4, 0))
        for label, w in [('#', 4), ('Position', 10),
                            ('Strand', 8), ('Guide (5\'3\')', 26),
                            ('PAM', 8), ('GC%', 6), ('Score', 8)]:
            tk.Label(header, text=label, bg=biosim_theme('panel'), fg='#0550ae',
                      font=('Arial', 10, 'bold'), width=w, anchor='w'
                      ).pack(side='left', padx=2, pady=4)

        # Rows
        for i, g in enumerate(guides, 1):
            color = ('#1a7f37' if g['score'] >= 80
                      else '#9a6700' if g['score'] >= 60
                      else '#bc4c00' if g['score'] >= 40
                      else '#cf222e')
            row = tk.Frame(self.crispr_results_inner, bg=biosim_theme('bg'))
            row.pack(fill='x', padx=10, pady=1)
            tk.Label(row, text=str(i), bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                      font=('Courier', 9), width=4, anchor='w'
                      ).pack(side='left', padx=2)
            tk.Label(row, text=str(g['position'] + 1),
                      bg=biosim_theme('bg'), fg=biosim_theme('text'),
                      font=('Courier', 9), width=10, anchor='w'
                      ).pack(side='left', padx=2)
            tk.Label(row, text=g['strand'],
                      bg=biosim_theme('bg'), fg='#8250df',
                      font=('Courier', 9, 'bold'), width=8, anchor='w'
                      ).pack(side='left', padx=2)
            tk.Label(row, text=g['guide'],
                      bg=biosim_theme('bg'), fg='#1a7f37',
                      font=('Courier', 10, 'bold'), width=26, anchor='w'
                      ).pack(side='left', padx=2)
            tk.Label(row, text=g['pam'],
                      bg=biosim_theme('bg'), fg='#9a6700',
                      font=('Courier', 9), width=8, anchor='w'
                      ).pack(side='left', padx=2)
            tk.Label(row, text='{:.0f}'.format(g['gc']),
                      bg=biosim_theme('bg'), fg=biosim_theme('text'),
                      font=('Courier', 9), width=6, anchor='w'
                      ).pack(side='left', padx=2)
            tk.Label(row, text='{:.0f}/100'.format(g['score']),
                      bg=biosim_theme('bg'), fg=color,
                      font=('Courier', 10, 'bold'), width=8, anchor='w'
                      ).pack(side='left', padx=2)
            # Copy button
            def make_copy(s):
                def copy():
                    self.clipboard_clear()
                    self.clipboard_append(s)
                    self.status_lbl.config(text='Guide copied',
                                                fg='#1a7f37')
                return copy
            tk.Button(row, text='Копировать', bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'),
                       relief='solid', bd=1, font=('Arial', 9), cursor='hand2',
                       padx=4, command=make_copy(g['guide'])
                       ).pack(side='left', padx=2)

        # Info footer
        tk.Label(self.crispr_results_inner,
                  text=('\nКритерии скоринга:\n'
                        '  • GC 40-60% (баланс эффективности)\n'
                        '  • G/C на 3\' конце (+10)\n'
                        '  • Нет TTTT (хорошая экспрессия с U6/T7)\n'
                        '  • Нет самокомплементарных участков\n\n'
                        'Это первичный дизайн. Проверь off-target в CRISPR-tools (e.g. CHOPCHOP, CRISPOR)'
                        if self.lang == 'ru' else
                        '\nScoring criteria:\n'
                        '  • GC 40-60% (efficiency balance)\n'
                        '  • G/C at 3\' end (+10)\n'
                        '  • No TTTT (good U6/T7 expression)\n'
                        '  • No self-complementary regions\n\n'
                        'Initial design only. Verify off-targets via CHOPCHOP, CRISPOR, etc.'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9, 'italic'), justify='left',
                  wraplength=1000).pack(anchor='w', padx=20, pady=10)

    def _build_tab_phylo(self):
        """Phylogenetic analysis: MSA + UPGMA tree."""
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=self.T['tab_phylo'])
        self.tab_phylo = tab

        # ─── Top controls ───
        ctrl = tk.LabelFrame(tab,
                                text='  ' + ('Филогенетический анализ'
                                                  if self.lang == 'ru' else
                                                  'Phylogenetic Analysis') + '  ',
                                bg=biosim_theme('bg'), fg='#8250df',
                                font=('Arial', 11, 'bold'), relief='flat')
        ctrl.pack(fill='x', padx=4, pady=4)

        tk.Label(ctrl,
                  text=('Введи последовательности в FASTA-формате (минимум 3 для дерева):'
                        if self.lang == 'ru' else
                        'Enter sequences in FASTA format (min 3 for tree):'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9, 'italic')).pack(anchor='w', padx=4, pady=2)

        seq_frame = tk.Frame(ctrl, bg=biosim_theme('panel'))
        seq_frame.pack(fill='x', padx=8, pady=4)
        sb_phylo = tk.Scrollbar(seq_frame)
        sb_phylo.pack(side='right', fill='y')
        self.phylo_text = tk.Text(seq_frame, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                                       insertbackground='#1c2128',
                                       relief='flat', font=('Courier', 9),
                                       height=8, padx=8, pady=4,
                                       yscrollcommand=sb_phylo.set)
        self.phylo_text.pack(side='left', fill='x', expand=True)
        sb_phylo.config(command=self.phylo_text.yview)

        # Example data (insulin from 4 species)
        example = ('>Human\nMALWMRLLPLLALLALWGPDPAAA\n'
                    '>Mouse\nMALLVHFLPLLALLAFWEPKPTSA\n'
                    '>Cow\nMALWTRLRPLLALLALWPPPPARA\n'
                    '>Chicken\nMALWIRSLPLLALLLLSGPGTSYA')
        self.phylo_text.insert('1.0', example)

        # Mode + protein checkbox + button
        opts_row = tk.Frame(ctrl, bg=biosim_theme('bg'))
        opts_row.pack(fill='x', padx=8, pady=4)

        self.phylo_protein_var = tk.BooleanVar(value=True)
        tk.Checkbutton(opts_row, text=('Белок (BLOSUM62)' if self.lang == 'ru'
                                          else 'Protein (BLOSUM62)'),
                        variable=self.phylo_protein_var,
                        bg=biosim_theme('bg'), fg='#0550ae',
                        selectcolor='#eef1f5',
                        activebackground='#f5f7fa',
                        font=('Arial', 9, 'bold')).pack(side='left', padx=4)

        tk.Button(opts_row, text='▶ ' + ('Построить дерево' if self.lang == 'ru'
                                              else 'Build tree'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2',
                   padx=14, pady=3, command=self._run_phylo
                   ).pack(side='right', padx=4)

        # ─── Results ───
        results_frame = tk.LabelFrame(tab,
                                          text='  ' + ('Дерево + MSA' if self.lang == 'ru'
                                                            else 'Tree + MSA') + '  ',
                                          bg=biosim_theme('bg'), fg='#1a7f37',
                                          font=('Arial', 11, 'bold'), relief='flat')
        results_frame.pack(fill='both', expand=True, padx=4, pady=4)

        canv = tk.Canvas(results_frame, bg=biosim_theme('bg'), highlightthickness=0)
        sb = tk.Scrollbar(results_frame, orient='vertical', command=canv.yview)
        canv.configure(yscrollcommand=sb.set)
        canv.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        sb.pack(side='right', fill='y')
        self.phylo_results_inner = tk.Frame(canv, bg=biosim_theme('bg'))
        self.phylo_results_inner.bind('<Configure>',
                                            lambda e: canv.configure(scrollregion=canv.bbox('all')))
        canv.create_window((0, 0), window=self.phylo_results_inner,
                             anchor='nw', width=1100)

    def _run_phylo(self):
        """Build phylogenetic tree from sequences."""
        # Parse FASTA
        text = self.phylo_text.get('1.0', 'end').strip()
        sequences = []
        labels = []
        current_label = None
        current_seq = []
        for line in text.split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current_label and current_seq:
                    sequences.append(''.join(current_seq).upper())
                    labels.append(current_label)
                current_label = line[1:].split()[0] if line[1:].strip() else 'unnamed'
                current_seq = []
            elif line:
                current_seq.append(''.join(c for c in line
                                             if c.isalpha() or c == '-'))
        if current_label and current_seq:
            sequences.append(''.join(current_seq).upper())
            labels.append(current_label)

        if len(sequences) < 3:
            for w in self.phylo_results_inner.winfo_children():
                w.destroy()
            tk.Label(self.phylo_results_inner,
                      text=('Нужно минимум 3 последовательности' if self.lang == 'ru'
                                      else 'Need at least 3 sequences'),
                      bg=biosim_theme('bg'), fg='#cf222e',
                      font=('Arial', 11, 'bold'),
                      padx=20, pady=20).pack(anchor='center')
            return

        is_prot = self.phylo_protein_var.get()

        # 1. Star alignment (MSA)
        msa = star_alignment(sequences, labels=labels, is_protein=is_prot)
        if msa is None:
            return
        aligned_seqs = msa['alignment']

        # 2. Distance matrix
        D, _ = distance_matrix_from_alignment(aligned_seqs, labels=labels)

        # 3. UPGMA tree
        tree = upgma_tree(D, labels)

        # Display
        for w in self.phylo_results_inner.winfo_children():
            w.destroy()

        # MSA section
        msa_frame = tk.LabelFrame(self.phylo_results_inner,
                                       text='  Multiple Sequence Alignment  ',
                                       bg=biosim_theme('bg'), fg='#0550ae',
                                       font=('Arial', 10, 'bold'), relief='flat')
        msa_frame.pack(fill='x', padx=10, pady=4)

        msa_text = tk.Text(msa_frame, bg=biosim_theme('panel'), fg=biosim_theme('text'),
                              font=('Courier', 10), relief='flat',
                              height=min(8, len(aligned_seqs) + 1),
                              wrap='none', padx=10, pady=6)
        msa_text.pack(fill='x', padx=4, pady=4)

        max_label = max(len(l) for l in msa['labels'])
        for label, seq in zip(msa['labels'], aligned_seqs):
            msa_text.insert('end', '{:<{w}}  {}\n'.format(label, seq[:80],
                                                              w=max_label))
        # Consensus
        cons = ''
        for i in range(len(aligned_seqs[0])):
            col = [s[i] for s in aligned_seqs if i < len(s)]
            non_gap = [c for c in col if c != '-']
            if non_gap and len(set(non_gap)) == 1:
                cons += non_gap[0]
            elif non_gap:
                cons += '.'
            else:
                cons += '-'
        msa_text.insert('end', '\n{:<{w}}  {}'.format('Consensus:', cons[:80],
                                                          w=max_label))
        msa_text.config(state='disabled')

        # Distance matrix
        dist_frame = tk.LabelFrame(self.phylo_results_inner,
                                        text='  Distance Matrix  ',
                                        bg=biosim_theme('bg'), fg='#9a6700',
                                        font=('Arial', 10, 'bold'), relief='flat')
        dist_frame.pack(fill='x', padx=10, pady=4)

        dist_inner = tk.Frame(dist_frame, bg=biosim_theme('bg'))
        dist_inner.pack(padx=10, pady=6)
        # Header
        tk.Label(dist_inner, text='', bg=biosim_theme('bg')).grid(row=0, column=0)
        for j, lbl in enumerate(labels):
            tk.Label(dist_inner, text=lbl[:8],
                      bg=biosim_theme('bg'), fg='#0550ae',
                      font=('Courier', 9, 'bold')
                      ).grid(row=0, column=j + 1, padx=4, pady=2)
        for i, lbl in enumerate(labels):
            tk.Label(dist_inner, text=lbl[:8],
                      bg=biosim_theme('bg'), fg='#0550ae',
                      font=('Courier', 9, 'bold'), anchor='e'
                      ).grid(row=i + 1, column=0, padx=4, pady=2, sticky='e')
            for j in range(len(labels)):
                val = D[i][j]
                if i == j:
                    color = '#57606a'
                elif val < 0.1:
                    color = '#1a7f37'
                elif val < 0.3:
                    color = '#9a6700'
                else:
                    color = '#cf222e'
                tk.Label(dist_inner, text='{:.3f}'.format(val),
                          bg=biosim_theme('bg'), fg=color,
                          font=('Courier', 9)
                          ).grid(row=i + 1, column=j + 1, padx=4, pady=2)

        # Tree visualization (ASCII)
        tree_frame = tk.LabelFrame(self.phylo_results_inner,
                                        text='  UPGMA Tree  ',
                                        bg=biosim_theme('bg'), fg='#8250df',
                                        font=('Arial', 10, 'bold'), relief='flat')
        tree_frame.pack(fill='x', padx=10, pady=4)

        # Visual tree with canvas
        tree_canvas = tk.Canvas(tree_frame, bg=biosim_theme('panel'), height=max(200, 30 * len(labels)),
                                    highlightthickness=0)
        tree_canvas.pack(fill='x', padx=4, pady=6)

        self._draw_phylo_tree(tree_canvas, tree, labels, D)

        # Newick text
        newick_frame = tk.Frame(self.phylo_results_inner, bg=biosim_theme('bg'))
        newick_frame.pack(fill='x', padx=10, pady=4)
        tk.Label(newick_frame, text='Newick format:',
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'),
                  font=('Arial', 9, 'bold')).pack(anchor='w', padx=4)
        newick_lbl = tk.Text(newick_frame, bg=biosim_theme('panel'), fg='#1a7f37',
                                 font=('Courier', 9), relief='flat',
                                 wrap='word', height=3, padx=10, pady=6)
        newick_lbl.pack(fill='x', padx=4, pady=2)
        newick_lbl.insert('1.0', tree['newick'])
        newick_lbl.config(state='disabled')

        # Copy button
        def copy_newick():
            self.clipboard_clear()
            self.clipboard_append(tree['newick'])
            self.status_lbl.config(text='Newick copied', fg='#1a7f37')
        tk.Button(newick_frame, text='Copy Newick',
                   bg=biosim_theme('panel'), highlightbackground='#1c2128', fg=biosim_theme('text'), relief='solid', bd=1,
                   font=('Arial', 9, 'bold'), cursor='hand2',
                   padx=10, pady=2, command=copy_newick
                   ).pack(anchor='e', padx=4, pady=2)

    def _draw_phylo_tree(self, canvas, tree, labels, D):
        """Draw phylogenetic tree on canvas (rectangular layout)."""
        canvas.delete('all')
        canvas.update_idletasks()
        w = canvas.winfo_width() if canvas.winfo_width() > 100 else 1000
        h = canvas.winfo_height() if canvas.winfo_height() > 100 else 200

        n = len(labels)
        if n < 2:
            return

        # Get max distance (depth)
        max_d = max(max(row) for row in D)
        if max_d <= 0:
            max_d = 1

        # Y positions for leaves
        y_step = (h - 40) / max(1, n - 1)
        leaf_positions = {}
        for i, label in enumerate(labels):
            leaf_positions[label] = 20 + i * y_step

        # X scale
        x_left = 20
        x_right = w - 200
        x_scale = (x_right - x_left) / max_d

        # Recursive walk merges to draw branches
        # Use merges from tree, in order — each merge joins 2 clusters at given distance
        cluster_pos = dict(leaf_positions)  # label -> y position

        # Walk merges
        for (a, b, dist) in tree['merges']:
            # a and b can be either leaf label or merged cluster (we don't have unique IDs)
            # Simplified: average y of all leaves in subtree
            def collect_leaves(label, originals):
                if label in originals:
                    return [label]
                # If label is a merged subtree, parse the Newick-like
                # For now, just match any leaf containing it
                found = [l for l in originals if l in label]
                return found if found else []

            leaves_a = collect_leaves(a, labels)
            leaves_b = collect_leaves(b, labels)
            if not leaves_a or not leaves_b:
                continue
            ya = sum(cluster_pos[l] for l in leaves_a if l in cluster_pos) / len(leaves_a)
            yb = sum(cluster_pos[l] for l in leaves_b if l in cluster_pos) / len(leaves_b)
            x_merge = x_left + (dist / 2) * x_scale

            # Vertical line connecting clusters
            canvas.create_line(x_merge, ya, x_merge, yb,
                                  fill='#1a7f37', width=2)
            # Horizontal lines from each cluster to merge point
            xa = x_left + (cluster_pos.get(a + '_x', 0)) if a in cluster_pos else x_left
            xb = x_left + (cluster_pos.get(b + '_x', 0)) if b in cluster_pos else x_left
            canvas.create_line(xa, ya, x_merge, ya, fill='#1a7f37', width=2)
            canvas.create_line(xb, yb, x_merge, yb, fill='#1a7f37', width=2)

            # Update cluster position
            new_label = '(' + a + ',' + b + ')'
            cluster_pos[new_label] = (ya + yb) / 2

        # Draw leaf labels and small branches to leaf
        for label, y in leaf_positions.items():
            canvas.create_line(x_left, y, x_left + 5, y, fill='#1a7f37', width=2)
            canvas.create_text(x_right + 10, y, text=label, anchor='w',
                                  fill='#0550ae', font=('Arial', 10, 'bold'))
            # Connect leaf x to some x_merge — draw horizontal line
            canvas.create_line(x_left, y, x_right, y, fill='#d0d7de',
                                  width=1, dash=(2, 4))

        # Scale bar
        canvas.create_text(x_left, h - 12, text='Distance ',
                              fill='#57606a', font=('Arial', 9, 'italic'),
                              anchor='w')
        canvas.create_line(x_left, h - 20, x_right, h - 20,
                              fill='#57606a', width=1, arrow='last')
        canvas.create_text(x_right, h - 12, text='{:.3f}'.format(max_d),
                              fill='#57606a', font=('Courier', 9), anchor='w')

    def _build_tab_msa(self):
        """Multiple sequence alignment tab."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text='MSA')
        tk.Label(tab, text=('Множественное выравнивание последовательностей'
                              if is_ru else 'Multiple Sequence Alignment'),
                  bg=biosim_theme('bg'), fg='#1a7f37', font=('Arial', 13, 'bold')).pack(pady=8)
        tk.Label(tab, text=('Вставьте последовательности — по одной на строку'
                              if is_ru else 'Paste sequences — one per line'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'), font=('Arial', 9)).pack()

        self.msa_input = scrolledtext.ScrolledText(
            tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'), font=('Courier', 10),
            height=6, relief='flat')
        self.msa_input.pack(fill='x', padx=16, pady=8)
        self.msa_input.insert('1.0',
            'ATGGCAAGCTTAGCTAGCTAGCATCG\nATGGCAAGCTAAGCTAGCTAGCATCG\n'
            'ATGGCTAGCTTAGCTAGGTAGCATCG\nATGGCAAGCTTAGCTACCTAGCATCG')

        self.msa_output = scrolledtext.ScrolledText(
            tab, bg=biosim_theme('bg'), fg=biosim_theme('text'), font=('Courier', 10),
            height=14, relief='flat', state='disabled')
        self.msa_output.pack(fill='both', expand=True, padx=16, pady=8)

        tk.Button(tab, text=('▶ Выровнять' if is_ru else '▶ Align'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=self._run_msa).pack(fill='x', padx=16, pady=(0, 8))

    def _run_msa(self):
        is_ru = (self.lang == 'ru')
        text = self.msa_input.get('1.0', 'end').strip()
        seqs = [s.strip() for s in text.split('\n') if s.strip()]
        if len(seqs) < 2:
            return
        try:
            r = msa_progressive(seqs)
            cons = msa_column_conservation(r['aligned'])
        except Exception as e:
            biosim_log('warn', '_run_msa: {}'.format(e))
            return
        out = self.msa_output
        out.config(state='normal'); out.delete('1.0', 'end')
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
        out.insert('end', '\n\n' + ('Идентичность: ' if is_ru else 'Identity: ')
                    + '{}%\n'.format(r['identity']), 'cons')
        out.config(state='disabled')

    def _build_tab_docking(self):
        """Molecular docking tab."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=('Докинг' if is_ru else 'Docking'))
        tk.Label(tab, text=('Молекулярный докинг (связывание с белком)'
                              if is_ru else 'Molecular Docking'),
                  bg=biosim_theme('bg'), fg='#0969da', font=('Arial', 13, 'bold')).pack(pady=8)

        row = tk.Frame(tab, bg=biosim_theme('bg')); row.pack(fill='x', padx=16, pady=8)
        tk.Label(row, text='SMILES:', bg=biosim_theme('bg'), fg=biosim_theme('text'),
                  font=('Arial', 9)).pack(side='left')
        self.dock_smiles_var = tk.StringVar(value='CC(=O)Oc1ccccc1C(=O)O')
        tk.Entry(row, textvariable=self.dock_smiles_var, bg=biosim_theme('panel2'),
                  fg=biosim_theme('text'), relief='flat').pack(side='left', fill='x',
                                                     expand=True, padx=6)

        self.dock_output = scrolledtext.ScrolledText(
            tab, bg=biosim_theme('bg'), fg=biosim_theme('text'), font=('Courier', 10),
            height=18, relief='flat', state='disabled')
        self.dock_output.pack(fill='both', expand=True, padx=16, pady=8)

        tk.Button(tab, text=('▶ Докинг' if is_ru else '▶ Dock'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=self._run_docking).pack(fill='x', padx=16, pady=(0, 8))

    def _run_docking(self):
        is_ru = (self.lang == 'ru')
        out = self.dock_output
        out.config(state='normal'); out.delete('1.0', 'end')
        for tag, c in [('ok', '#1a7f37'), ('bad', '#cf222e'), ('hdr', '#0969da')]:
            out.tag_config(tag, foreground=c)
        if not _ensure_rdkit_p4():
            out.insert('end', ('RDKit не установлен\n' if is_ru else 'RDKit not installed\n'), 'bad')
            out.config(state='disabled'); return
        out.insert('end', ('⏳ Генерация 3D и докинг...\n' if is_ru else '⏳ Docking...\n'))
        out.config(state='disabled'); out.update()
        try:
            r = dock_smiles(self.dock_smiles_var.get().strip())
        except Exception as e:
            biosim_log('warn', '_run_docking: {}'.format(e))
            r = {'error': str(e)[:80]}
        out.config(state='normal'); out.delete('1.0', 'end')
        for tag, c in [('ok', '#1a7f37'), ('bad', '#cf222e'), ('hdr', '#0969da')]:
            out.tag_config(tag, foreground=c)
        if 'error' in r:
            out.insert('end', '' + r['error'] + '\n', 'bad')
            out.config(state='disabled'); return
        out.insert('end', ('Результат докинга:\n' if is_ru else 'Docking result:\n'), 'hdr')
        out.insert('end', '  ' + ('Атомов лиганда: ' if is_ru else 'Ligand atoms: ')
                    + str(r['ligand_atoms']) + '\n')
        out.insert('end', '  ' + ('Проверено поз: ' if is_ru else 'Poses: ')
                    + str(r['poses_tested']) + '\n\n')
        score = r['best_score']
        out.insert('end', '  ' + ('Лучший score: ' if is_ru else 'Best score: ')
                    + str(score) + '\n', 'ok' if score < 0 else 'bad')
        d = r.get('details', {})
        if d:
            out.insert('end', '  ' + ('Контактов: ' if is_ru else 'Contacts: ')
                        + str(d.get('contacts', 0)) + '\n')
            out.insert('end', '  ' + ('Водородных связей: ' if is_ru else 'H-bonds: ')
                        + str(d.get('h_bonds', 0)) + '\n')
            out.insert('end', '  ' + ('Качество: ' if is_ru else 'Quality: ')
                        + str(d.get('binding_quality', '?')) + '\n', 'ok')
        kd = r.get('estimated_kd_nm', 0)
        if kd < 99999:
            out.insert('end', '\n  ' + ('Прибл. Kd: ' if is_ru else 'Est. Kd: ')
                        + '{} nM\n'.format(kd))
        out.config(state='disabled')

    def _build_tab_adme(self):
        """ADME / drug-likeness tab."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text='ADME')
        tk.Label(tab, text=('Анализ лекарственности (ADME)'
                              if is_ru else 'Drug-likeness (ADME)'),
                  bg=biosim_theme('bg'), fg='#8250df', font=('Arial', 13, 'bold')).pack(pady=8)

        row = tk.Frame(tab, bg=biosim_theme('bg')); row.pack(fill='x', padx=16, pady=8)
        tk.Label(row, text='SMILES:', bg=biosim_theme('bg'), fg=biosim_theme('text'),
                  font=('Arial', 9)).pack(side='left')
        self.adme_smiles_var = tk.StringVar(value='CC(=O)Oc1ccccc1C(=O)O')
        tk.Entry(row, textvariable=self.adme_smiles_var, bg=biosim_theme('panel2'),
                  fg=biosim_theme('text'), relief='flat').pack(side='left', fill='x',
                                                     expand=True, padx=6)

        self.adme_output = scrolledtext.ScrolledText(
            tab, bg=biosim_theme('bg'), fg=biosim_theme('text'), font=('Arial', 10),
            height=18, relief='flat', state='disabled')
        self.adme_output.pack(fill='both', expand=True, padx=16, pady=8)

        tk.Button(tab, text=('▶ Анализировать' if is_ru else '▶ Analyze'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat',
                   font=('Arial', 10, 'bold'), cursor='hand2', pady=4,
                   command=self._run_adme_tab).pack(fill='x', padx=16, pady=(0, 8))

    def _run_adme_tab(self):
        is_ru = (self.lang == 'ru')
        out = self.adme_output
        out.config(state='normal'); out.delete('1.0', 'end')
        for tag, c in [('ok', '#1a7f37'), ('bad', '#cf222e'),
                        ('hdr', '#0969da'), ('dim', '#57606a')]:
            out.tag_config(tag, foreground=c,
                            font=('Arial', 10, 'bold') if tag == 'hdr' else ('Arial', 10))
        if not _ensure_rdkit_p4():
            out.insert('end', ('RDKit не установлен\n' if is_ru else 'RDKit not installed\n'), 'bad')
            out.config(state='disabled'); return
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski
            mol = Chem.MolFromSmiles(self.adme_smiles_var.get().strip())
            if mol is None:
                out.insert('end', ('Неверный SMILES\n' if is_ru else 'Invalid SMILES\n'), 'bad')
                out.config(state='disabled'); return
            mw = Descriptors.MolWt(mol); logp = Descriptors.MolLogP(mol)
            hd = Lipinski.NumHDonors(mol); ha = Lipinski.NumHAcceptors(mol)
            rb = Lipinski.NumRotatableBonds(mol); tpsa = Descriptors.TPSA(mol)
            r = adme_properties(mw, logp, hd, ha, rb, tpsa)
        except Exception as e:
            biosim_log('warn', '_run_adme_tab: {}'.format(e))
            out.insert('end', '' + str(e)[:80] + '\n', 'bad')
            out.config(state='disabled'); return
        out.insert('end', ('Свойства молекулы:\n' if is_ru else 'Properties:\n'), 'hdr')
        out.insert('end', '  MW = {:.1f} | LogP = {:.2f}\n'.format(mw, logp))
        out.insert('end', '  H-доноры={} H-акцепторы={}\n'.format(hd, ha))
        out.insert('end', '  Вращ.связи={} TPSA={:.1f}\n\n'.format(rb, tpsa))
        out.insert('end', ('Правило Липински:\n' if is_ru else 'Lipinski:\n'), 'hdr')
        for rule, ok in r['lipinski'].items():
            out.insert('end', '  {} {}\n'.format('' if ok else '', rule),
                        'ok' if ok else 'bad')
        out.insert('end', ('\nПравило Вебера:\n' if is_ru else '\nVeber:\n'), 'hdr')
        for rule, ok in r['veber'].items():
            out.insert('end', '  {} {}\n'.format('' if ok else '', rule),
                        'ok' if ok else 'bad')
        out.insert('end', '\n  {} {}\n'.format('' if r['druglike'] else '',
                    ('Лекарственноподобна' if is_ru else 'Drug-like')),
                    'ok' if r['druglike'] else 'bad')
        out.insert('end', '  {} {}\n'.format('' if r['oral_bioavailable'] else '',
                    ('Перорально биодоступна' if is_ru else 'Orally bioavailable')),
                    'ok' if r['oral_bioavailable'] else 'bad')
        out.insert('end', '  {} {}\n'.format('' if r['bbb_permeable'] else '',
                    ('Проникает в мозг' if is_ru else 'BBB permeable')),
                    'ok' if r['bbb_permeable'] else 'dim')
        out.config(state='disabled')


    def _build_tab_orf(self):
        """Advanced ORF finder tab."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text='ORF')
        tk.Label(tab, text=('Поиск открытых рамок считывания (6 рамок)'
                              if is_ru else 'ORF Finder (6 frames)'),
                  bg=biosim_theme('bg'), fg='#1a7f37', font=('Arial', 13, 'bold')).pack(pady=8)
        tk.Label(tab, text=('Вставьте ДНК последовательность'
                              if is_ru else 'Paste DNA sequence'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'), font=('Arial', 9)).pack()
        self.orf_input = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                     font=('Courier', 10), height=5, relief='flat')
        self.orf_input.pack(fill='x', padx=16, pady=8)
        self.orf_input.insert('1.0', 'ATGAAACTGATCCTGGCAGCATAAGGGCCCATGGCATTTTGA')
        self.orf_output = scrolledtext.ScrolledText(tab, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                                     font=('Courier', 10), height=15, relief='flat', state='disabled')
        self.orf_output.pack(fill='both', expand=True, padx=16, pady=8)
        tk.Button(tab, text=('▶ Найти ORF' if is_ru else '▶ Find ORFs'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=self._run_orf_finder).pack(fill='x', padx=16, pady=(0, 8))

    def _run_orf_finder(self):
        is_ru = (self.lang == 'ru')
        dna = self.orf_input.get('1.0', 'end').strip()
        orfs = find_orfs_advanced(dna, min_length=15)
        out = self.orf_output
        out.config(state='normal'); out.delete('1.0', 'end')
        out.tag_config('hdr', foreground='#0969da', font=('Courier', 10, 'bold'))
        out.tag_config('orf', foreground='#1a7f37')
        out.insert('end', ('Найдено ORF: ' if is_ru else 'ORFs found: ') + str(len(orfs)) + '\n\n', 'hdr')
        for i, orf in enumerate(orfs[:20], 1):
            out.insert('end', '#{} '.format(i), 'hdr')
            out.insert('end', '{} {} {}-{} ({} аа)\n'.format(
                ('цепь' if is_ru else 'strand'), orf['strand'],
                orf['start'], orf['end'], orf['length_aa']))
            out.insert('end', '   ' + orf['protein'][:50] + ('...' if len(orf['protein']) > 50 else '') + '\n', 'orf')
        out.config(state='disabled')

    def _build_tab_rna(self):
        """RNA secondary structure tab."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text='RNA')
        tk.Label(tab, text=('Вторичная структура РНК (Nussinov)'
                              if is_ru else 'RNA secondary structure (Nussinov)'),
                  bg=biosim_theme('bg'), fg='#0969da', font=('Arial', 13, 'bold')).pack(pady=8)
        tk.Label(tab, text=('Вставьте РНК/ДНК (макс 300 нт)'
                              if is_ru else 'Paste RNA/DNA (max 300 nt)'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'), font=('Arial', 9)).pack()
        self.rna_input = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                    font=('Courier', 10), height=4, relief='flat')
        self.rna_input.pack(fill='x', padx=16, pady=8)
        self.rna_input.insert('1.0', 'GGGAAAUUUCCCGGGAAACCC')
        # Graphical canvas for arc diagram (like professional RNA tools)
        self.rna_canvas = tk.Canvas(tab, bg=biosim_theme('panel'), height=170,
                                     relief='flat', highlightthickness=1,
                                     highlightbackground=biosim_theme('border'))
        self.rna_canvas.pack(fill='x', padx=16, pady=(0, 6))
        self.rna_output = scrolledtext.ScrolledText(tab, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                                    font=('Courier', 11), height=8, relief='flat', state='disabled')
        self.rna_output.pack(fill='both', expand=True, padx=16, pady=8)
        tk.Button(tab, text=('▶ Свернуть' if is_ru else '▶ Fold'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=self._run_rna_fold).pack(fill='x', padx=16, pady=(0, 8))

    def _draw_rna_arcs(self, seq, pairs):
        """Draw an arc diagram of RNA base pairs (professional visualization).
        Each base is a dot on a line; pairs are connected by coloured arcs."""
        cv = self.rna_canvas
        cv.delete('all')
        n = len(seq)
        if n == 0:
            return
        cv.update_idletasks()
        W = cv.winfo_width() or 600
        margin = 24
        usable = max(50, W - 2 * margin)
        baseline = 150
        step = usable / max(1, n - 1) if n > 1 else usable
        # Colour per nucleotide
        nt_color = {'A': '#1a7f37', 'U': '#cf222e', 'T': '#cf222e',
                    'G': '#0969da', 'C': '#9a6700'}
        # Draw arcs first (so dots sit on top)
        paired_set = set()
        for (i, j) in pairs:
            paired_set.add(i); paired_set.add(j)
            x1 = margin + i * step
            x2 = margin + j * step
            r = (x2 - x1) / 2
            cx = (x1 + x2) / 2
            # Arc height proportional to span
            h = min(120, 20 + r * 0.7)
            cv.create_arc(cx - r, baseline - h, cx + r, baseline + h,
                          start=0, extent=180, style='arc',
                          outline='#8250df', width=2)
        # Draw baseline
        cv.create_line(margin, baseline, margin + usable, baseline,
                       fill='#d0d7de', width=1)
        # Draw nucleotide dots + letters
        for k, base in enumerate(seq):
            x = margin + k * step
            col = nt_color.get(base, '#57606a')
            # Highlight paired bases with filled dot
            if k in paired_set:
                cv.create_oval(x - 4, baseline - 4, x + 4, baseline + 4,
                               fill=col, outline=col)
            else:
                cv.create_oval(x - 3, baseline - 3, x + 3, baseline + 3,
                               fill='#ffffff', outline=col, width=1)
            # Letter labels only if not too crowded
            if step > 9:
                cv.create_text(x, baseline + 14, text=base, fill=col,
                               font=('Courier', 8))

    def _run_rna_fold(self):
        is_ru = (self.lang == 'ru')
        rna = self.rna_input.get('1.0', 'end').strip()
        r = rna_secondary_structure(rna)
        out = self.rna_output
        out.config(state='normal'); out.delete('1.0', 'end')
        out.tag_config('hdr', foreground='#0969da', font=('Courier', 11, 'bold'))
        out.tag_config('seq', foreground=biosim_theme('text'), font=('Courier', 11))
        out.tag_config('str', foreground='#1a7f37', font=('Courier', 11))
        out.tag_config('ok', foreground='#1a7f37', font=('Courier', 10, 'bold'))
        if r.get('error'):
            out.insert('end', '' + r['error'] + '\n', 'hdr')
            out.config(state='disabled')
            self.rna_canvas.delete('all')
            return
        seq = r['sequence']; struct = r['structure']
        pairs = r.get('pairs', [])
        # Draw the graphical arc diagram
        try:
            self._draw_rna_arcs(seq, pairs)
        except Exception as _exc:
            biosim_log('warn', '_run_rna_fold: {}'.format(_exc))
        # Text output below
        out.insert('end', ('✓ Структура предсказана\n\n' if is_ru
                            else '✓ Structure predicted\n\n'), 'ok')
        out.insert('end', ('Последовательность и структура:\n\n' if is_ru else 'Sequence and structure:\n\n'), 'hdr')
        # Show in chunks of 50
        for i in range(0, len(seq), 50):
            out.insert('end', seq[i:i+50] + '\n', 'seq')
            out.insert('end', struct[i:i+50] + '\n\n', 'str')
        npairs = r.get('num_pairs', len(pairs))
        out.insert('end', ('Пар оснований: ' if is_ru else 'Base pairs: ') + str(npairs) + '\n', 'hdr')
        out.insert('end', '( = ' + ('спаренное' if is_ru else 'paired') + ', . = ' + ('свободное' if is_ru else 'unpaired') + '\n', 'seq')
        out.config(state='disabled')

    def _build_tab_protein(self):
        """Protein properties tab."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=('Белок' if is_ru else 'Protein'))
        tk.Label(tab, text=('Физико-химические свойства белка'
                              if is_ru else 'Protein physicochemical properties'),
                  bg=biosim_theme('bg'), fg='#8250df', font=('Arial', 13, 'bold')).pack(pady=8)
        tk.Label(tab, text=('Вставьте белок (1-буквенный код)'
                              if is_ru else 'Paste protein (1-letter)'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'), font=('Arial', 9)).pack()
        self.prot_input = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                     font=('Courier', 10), height=5, relief='flat')
        self.prot_input.pack(fill='x', padx=16, pady=8)
        self.prot_input.insert('1.0', 'MAKVLILACGDEFGHIKLMNPQRSTVWY')
        # Graphical hydropathy plot frame (Kyte-Doolittle, like Vector NTI/ExPASy)
        self.prot_fig_frame = tk.Frame(tab, bg=biosim_theme('bg'))
        self.prot_fig_frame.pack(fill='x', padx=16, pady=(0, 4))
        self._prot_canvas = {'c': None}
        self.prot_output = scrolledtext.ScrolledText(tab, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                                     font=('Arial', 10), height=10, relief='flat', state='disabled')
        self.prot_output.pack(fill='both', expand=True, padx=16, pady=8)
        tk.Button(tab, text=('▶ Анализировать' if is_ru else '▶ Analyze'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=self._run_protein_props).pack(fill='x', padx=16, pady=(0, 8))

    def _run_protein_props(self):
        is_ru = (self.lang == 'ru')
        protein = self.prot_input.get('1.0', 'end').strip()
        r = protein_properties(protein)
        out = self.prot_output
        out.config(state='normal'); out.delete('1.0', 'end')
        out.tag_config('hdr', foreground='#0969da', font=('Arial', 11, 'bold'))
        out.tag_config('val', foreground='#1a7f37', font=('Arial', 10, 'bold'))
        out.tag_config('ok', foreground='#1a7f37', font=('Arial', 10, 'bold'))
        if r.get('error'):
            out.insert('end', '' + r['error'] + '\n'); out.config(state='disabled'); return
        # Draw Kyte-Doolittle hydropathy plot
        try:
            self._draw_hydropathy(protein.upper())
        except Exception as _exc:
            biosim_log('warn', '_run_protein_props: {}'.format(_exc))
        out.insert('end', ('✓ Анализ завершён\n\n' if is_ru else '✓ Analysis complete\n\n'), 'ok')
        out.insert('end', ('Свойства белка:\n\n' if is_ru else 'Protein properties:\n\n'), 'hdr')
        rows = [
            (('Длина' if is_ru else 'Length'), '{} аа'.format(r['length'])),
            (('Молек. масса' if is_ru else 'MW'), '{} kDa ({} Da)'.format(r['mw_kda'], r['mw_da'])),
            (('Изоэлектр. точка pI' if is_ru else 'Isoelectric pI'), str(r['pI'])),
            ('GRAVY', '{} ({})'.format(r['gravy'], ('гидрофобный' if r['hydrophobic'] else 'гидрофильный') if is_ru else ('hydrophobic' if r['hydrophobic'] else 'hydrophilic'))),
            (('Положит. заряд' if is_ru else 'Positive charge'), '+{}'.format(r['charged_pos'])),
            (('Отрицат. заряд' if is_ru else 'Negative charge'), '-{}'.format(r['charged_neg'])),
            (('Чистый заряд' if is_ru else 'Net charge'), str(r['net_charge'])),
        ]
        for label, val in rows:
            out.insert('end', '  {}: '.format(label))
            out.insert('end', val + '\n', 'val')
        out.config(state='disabled')

    def _draw_hydropathy(self, protein):
        """Kyte-Doolittle hydropathy plot — classic protein hydrophobicity
        profile (used in Vector NTI / ExPASy ProtScale). Peaks above 0 =
        hydrophobic regions (possible membrane-spanning)."""
        kd = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
              'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
              'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
              'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
        seq = [a for a in protein if a in kd]
        if len(seq) < 5:
            return
        window = min(9, len(seq))
        half = window // 2
        xs, ys = [], []
        for i in range(len(seq)):
            lo = max(0, i - half); hi = min(len(seq), i + half + 1)
            vals = [kd[a] for a in seq[lo:hi]]
            xs.append(i + 1); ys.append(sum(vals) / len(vals))
        is_ru = (self.lang == 'ru')
        try:
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if self._prot_canvas['c']:
                self._prot_canvas['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(7, 2.6), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            ax.plot(xs, ys, color='#8250df', linewidth=1.5)
            ax.axhline(y=0, color='#57606a', linestyle='-', alpha=0.4)
            ax.fill_between(xs, ys, 0, where=[y > 0 for y in ys],
                            color='#cf222e', alpha=0.25)
            ax.fill_between(xs, ys, 0, where=[y <= 0 for y in ys],
                            color='#0969da', alpha=0.25)
            ax.set_xlabel(('Позиция (аа)' if is_ru else 'Position (aa)'), color='#1c2128')
            ax.set_ylabel(('Гидрофобность' if is_ru else 'Hydropathy'), color='#1c2128')
            ax.set_title(('Профиль Kyte-Doolittle' if is_ru else 'Kyte-Doolittle profile'),
                         color='#1c2128', fontsize=10)
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            fig.tight_layout()
            cv = FigureCanvasTkAgg(fig, master=self.prot_fig_frame); cv.draw()
            cv.get_tk_widget().pack(fill='x'); self._prot_canvas['c'] = cv
        except Exception as _exc:
            biosim_log('warn', '_draw_hydropathy: {}'.format(_exc))

    def _build_tab_gcanalysis(self):
        """GC content and skew analysis tab."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text='GC')
        tk.Label(tab, text=('GC-состав, GC-skew, CAI'
                              if is_ru else 'GC content, GC skew, CAI'),
                  bg=biosim_theme('bg'), fg='#9a6700', font=('Arial', 13, 'bold')).pack(pady=8)
        self.gc_input = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                   font=('Courier', 10), height=5, relief='flat')
        self.gc_input.pack(fill='x', padx=16, pady=8)
        self.gc_input.insert('1.0', 'ATGAAACTGATCCTGGCAGCAGGGCCCGCGCGCGCATATATATGCGCGCATGCAT')
        self.gc_fig_frame = tk.Frame(tab, bg=biosim_theme('bg'))
        self.gc_fig_frame.pack(fill='both', expand=True, padx=16, pady=8)
        self.gc_info = tk.Label(tab, text='', bg=biosim_theme('bg'), fg='#1a7f37', font=('Arial', 10, 'bold'))
        self.gc_info.pack()
        self._gc_canvas = {'c': None}
        tk.Button(tab, text=('▶ Анализ' if is_ru else '▶ Analyze'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=self._run_gc_analysis).pack(fill='x', padx=16, pady=(0, 8))

    def _run_gc_analysis(self):
        is_ru = (self.lang == 'ru')
        dna = self.gc_input.get('1.0', 'end').strip()
        if len(dna) < 20:
            self.gc_info.config(
                text=('⚠ Нужно ≥20 нуклеотидов' if is_ru else '⚠ Need ≥20 nt'),
                fg='#cf222e')
            return
        gc_win = gc_content_windowed(dna, window=min(20, len(dna) // 3))
        cai = calculate_cai(dna)
        # Compute cumulative GC skew (G-C)/(G+C) — finds replication origin
        seq = dna.upper()
        skew = []
        cum = 0
        for i, base in enumerate(seq):
            if base == 'G':
                cum += 1
            elif base == 'C':
                cum -= 1
            skew.append(cum)
        try:
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            if self._gc_canvas['c']:
                self._gc_canvas['c'].get_tk_widget().destroy()
            fig = plt.Figure(figsize=(7, 3.6), facecolor='#f5f7fa')
            ax = fig.add_subplot(111, facecolor='#ffffff')
            if gc_win:
                xs = [g[0] for g in gc_win]; ys = [g[1] for g in gc_win]
                ax.plot(xs, ys, color='#9a6700', marker='o', markersize=3,
                        label='GC %' )
                ax.axhline(y=50, color='#57606a', linestyle='--', alpha=0.4)
            ax.set_xlabel(('Позиция' if is_ru else 'Position'), color='#1c2128')
            ax.set_ylabel('GC %', color='#9a6700')
            ax.tick_params(colors='#57606a'); ax.grid(alpha=0.2)
            for s in ax.spines.values(): s.set_color('#d0d7de')
            # Second axis: cumulative GC skew
            ax2 = ax.twinx()
            ax2.plot(range(len(skew)), skew, color='#0969da', alpha=0.7,
                     linewidth=1.5, label='GC skew')
            ax2.set_ylabel(('GC skew (накопл.)' if is_ru else 'GC skew (cumul.)'),
                           color='#0969da')
            ax2.tick_params(axis='y', colors='#0969da')
            fig.tight_layout()
            cv = FigureCanvasTkAgg(fig, master=self.gc_fig_frame); cv.draw()
            cv.get_tk_widget().pack(fill='both', expand=True); self._gc_canvas['c'] = cv
        except Exception as _exc:
            biosim_log('warn', '_run_gc_analysis: {}'.format(_exc))
        total_gc = (dna.upper().count('G') + dna.upper().count('C')) / len(dna) * 100
        self.gc_info.config(text='✓ GC = {:.1f}% | CAI = {} ({})'.format(
            total_gc, cai['cai'], cai.get('interpretation', '')), fg='#1a7f37')

    def _build_tab_digest(self):
        """Restriction digest tab."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=('Дайджест' if is_ru else 'Digest'))
        tk.Label(tab, text=('Виртуальный рестрикционный дайджест'
                              if is_ru else 'Virtual restriction digest'),
                  bg=biosim_theme('bg'), fg='#cf222e', font=('Arial', 13, 'bold')).pack(pady=8)
        self.digest_input = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                       font=('Courier', 10), height=5, relief='flat')
        self.digest_input.pack(fill='x', padx=16, pady=8)
        self.digest_input.insert('1.0', 'GAATTCAAAAAGGATCCTTTTTAAGCTTCCCCCGTCGACGGGGG')
        self.digest_output = scrolledtext.ScrolledText(tab, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                                       font=('Courier', 10), height=14, relief='flat', state='disabled')
        self.digest_output.pack(fill='both', expand=True, padx=16, pady=8)
        tk.Button(tab, text=('▶ Разрезать' if is_ru else '▶ Digest'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white', relief='flat', font=('Arial', 10, 'bold'),
                   cursor='hand2', pady=4, command=self._run_digest).pack(fill='x', padx=16, pady=(0, 8))

    def _build_tab_blast_pro(self):
        """BLAST search tab inside MolBio Workspace."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text='BLAST')
        tk.Label(tab, text=('BLAST — поиск похожих последовательностей (e-value, bit-score)'
                              if is_ru else 'BLAST — find similar sequences (e-value, bit-score)'),
                  bg=biosim_theme('bg'), fg='#0969da', font=('Arial', 12, 'bold')).pack(pady=6)
        tk.Label(tab, text=('Запрос:' if is_ru else 'Query:'), bg=biosim_theme('bg'),
                  fg=biosim_theme('text'), font=('Arial', 9)).pack(anchor='w', padx=14)
        self.blast_query = tk.Entry(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'), relief='flat', font=('Courier', 10))
        self.blast_query.pack(fill='x', padx=14, pady=3)
        self.blast_query.insert(0, 'ATGGCAAGCTTAGCTAGC')
        tk.Label(tab, text=('База (по одной на строку):' if is_ru else 'Database (one per line):'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text'), font=('Arial', 9)).pack(anchor='w', padx=14)
        self.blast_db = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                   font=('Courier', 9), height=4, relief='flat')
        self.blast_db.pack(fill='x', padx=14, pady=3)
        self.blast_db.insert('1.0', 'ATGGCAAGCTTAGCTAGCTAGCATCG\nTTTGGGCCCAAATTTGGG\nATGGCAAGCTTACCCGGGTAG')
        self.blast_results = scrolledtext.ScrolledText(tab, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                                       font=('Courier', 9), height=14, relief='flat', state='disabled')
        self.blast_results.pack(fill='both', expand=True, padx=14, pady=6)
        tk.Button(tab, text=('▶ BLAST поиск' if is_ru else '▶ BLAST search'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=self._run_blast_pro).pack(fill='x', padx=14, pady=4)

    def _run_blast_pro(self):
        is_ru = (self.lang == 'ru')
        query = self.blast_query.get().strip()
        db = [s.strip() for s in self.blast_db.get('1.0', 'end').split('\n') if s.strip()]
        if not query or not db:
            return
        try:
            hits = blast_search(query, db)
        except Exception as e:
            biosim_log('warn', '_run_blast_pro: {}'.format(e))
            return
        o = self.blast_results; o.config(state='normal'); o.delete('1.0', 'end')
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

    def _build_tab_annotation_pro(self):
        """Plasmid annotation tab inside MolBio Workspace."""
        is_ru = (self.lang == 'ru')
        tab = tk.Frame(self.notebook, bg=biosim_theme('bg'))
        self.notebook.add(tab, text=('Аннотация' if is_ru else 'Annotation'))
        tk.Label(tab, text=('Авто-аннотация плазмиды (промоторы, гены, ori, маркеры)'
                              if is_ru else 'Auto-annotate plasmid (promoters, genes, ori, markers)'),
                  bg=biosim_theme('bg'), fg='#1a7f37', font=('Arial', 12, 'bold')).pack(pady=6)
        self.annot_input = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                      font=('Courier', 9), height=5, relief='flat')
        self.annot_input.pack(fill='x', padx=14, pady=6)
        self.annot_input.insert('1.0', 'TAATACGACTCACTATAGGGAGGAGGATG' + 'GCAGCT' * 60 + 'TAA')
        self.annot_results = scrolledtext.ScrolledText(tab, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                                       font=('Courier', 9), height=15, relief='flat', state='disabled')
        self.annot_results.pack(fill='both', expand=True, padx=14, pady=6)
        tk.Button(tab, text=('▶ Аннотировать' if is_ru else '▶ Annotate'), bg='#1a7f37', highlightbackground='#1a7f37',
                   fg='white', relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=self._run_annotation_pro).pack(fill='x', padx=14, pady=4)

    def _run_annotation_pro(self):
        is_ru = (self.lang == 'ru')
        seq = self.annot_input.get('1.0', 'end').strip()
        if len(seq) < 30:
            return
        try:
            ann = annotate_plasmid(seq, min_orf_len=60)
        except Exception as _exc:
            biosim_log('warn', '_run_annotation_pro: {}'.format(_exc))
            return
        o = self.annot_results; o.config(state='normal'); o.delete('1.0', 'end')
        o.tag_config('hdr', foreground='#0969da', font=('Courier', 9, 'bold'))
        colors = {'promoter': '#9a6700', 'ORF': '#1a7f37', 'RBS': '#0550ae',
                  'ori': '#cf222e', 'resistance': '#8250df'}
        for t, c in colors.items():
            o.tag_config(t, foreground=c)
        o.insert('end', 'Length: {} bp | GC: {}% | {} features\n\n'.format(
            ann['length'], ann['gc_content'], ann['feature_count']), 'hdr')
        for f in ann['features']:
            o.insert('end', '{:>5}-{:<5} [{}] '.format(f['start'], f['end'], f['strand']))
            o.insert('end', '{}: {}\n'.format(f['type'], f['name']),
                      f['type'] if f['type'] in colors else 'hdr')
        o.config(state='disabled')


    def _run_digest(self):
        is_ru = (self.lang == 'ru')
        dna = self.digest_input.get('1.0', 'end').strip()
        r = restriction_digest(dna)
        out = self.digest_output
        out.config(state='normal'); out.delete('1.0', 'end')
        out.tag_config('hdr', foreground='#0969da', font=('Courier', 10, 'bold'))
        out.tag_config('cut', foreground='#cf222e')
        out.tag_config('frag', foreground='#1a7f37')
        out.insert('end', ('Разрезов: ' if is_ru else 'Cuts: ') + str(r['n_cuts']) + '\n', 'hdr')
        out.insert('end', ('Ферменты: ' if is_ru else 'Enzymes: ') + ', '.join(r['enzymes_that_cut']) + '\n\n', 'hdr')
        out.insert('end', ('Сайты разрезания:\n' if is_ru else 'Cut sites:\n'), 'hdr')
        for cut in r['cuts']:
            out.insert('end', '  {} @ {} ({})\n'.format(cut['enzyme'], cut['position'], cut['site']), 'cut')
        out.insert('end', '\n' + ('Фрагменты (bp):\n' if is_ru else 'Fragments (bp):\n'), 'hdr')
        for frag in r['fragments']:
            if frag > 0:
                bar = '█' * min(40, max(1, frag // 2))
                out.insert('end', '  {:>5} {}\n'.format(frag, bar), 'frag')
        out.config(state='disabled')


    def _struct_show_error(self, msg):
        self._struct_clear_results()
        tk.Label(self.struct_results_inner, text='' + msg,
                  bg=biosim_theme('bg'), fg='#cf222e',
                  font=('Arial', 11, 'bold'),
                  padx=20, pady=20).pack(anchor='center')


    # ─── Quick actions ───

    def _quick_paste(self):
        """Paste sequence from clipboard."""
        try:
            text = self.clipboard_get()
        except Exception as _exc:
            biosim_log('warn', '_quick_paste: {}'.format(_exc))
            self.status_lbl.config(text='Clipboard empty', fg='#cf222e')
            return
        if not text:
            return
        # Clean: remove whitespace, numbers, etc
        cleaned = ''.join(c.upper() for c in text if c.isalpha())
        self.current_sequence = cleaned
        self.sequence_type = detect_sequence_type(cleaned)
        self._update_status()
        # Reload into sequence editor if open
        if hasattr(self, 'seq_text'):
            try:
                self.seq_text.delete('1.0', 'end')
                self.seq_text.insert('1.0', cleaned)
                if hasattr(self, 'seq_type_var'):
                    self.seq_type_var.set(self.sequence_type)
                self._refresh_seq_analysis()
                self._refresh_seq_highlights()
                self.notebook.select(self.tab_sequence)
            except Exception as _exc:
                biosim_log('warn', '_quick_paste: {}'.format(_exc))
        self.status_lbl.config(text=('Вставлено · ' if self.lang == 'ru'
                                       else 'Pasted · ')
                                    + '{} bp · '.format(len(cleaned))
                                    + self.sequence_type.upper(),
                                fg='#1a7f37')

    def _quick_open_file(self):
        """Open FASTA or GenBank file."""
        from tkinter import filedialog
        filepath = filedialog.askopenfilename(
            title='Open sequence file',
            filetypes=[('FASTA', '*.fasta *.fa *.fna *.ffn *.faa'),
                        ('GenBank', '*.gb *.gbk *.genbank'),
                        ('Text', '*.txt'),
                        ('All files', '*.*')])
        if not filepath:
            return
        try:
            seq, name, seq_type = self._parse_sequence_file(filepath)
            if seq:
                self.current_sequence = seq
                self.sequence_name = name
                self.sequence_type = seq_type
                self._update_status()
                # Reload into editor
                if hasattr(self, 'seq_text'):
                    try:
                        self.seq_text.delete('1.0', 'end')
                        self.seq_text.insert('1.0', seq)
                        if hasattr(self, 'seq_name_var'):
                            self.seq_name_var.set(name)
                        if hasattr(self, 'seq_type_var'):
                            self.seq_type_var.set(seq_type)
                        self._refresh_seq_analysis()
                        self._refresh_seq_highlights()
                        self.notebook.select(self.tab_sequence)
                    except Exception as _exc:
                        biosim_log('warn', '_quick_open_file: {}'.format(_exc))
                self.status_lbl.config(
                    text='' + name + ' · ' + '{} bp'.format(len(seq)),
                    fg='#1a7f37')
            else:
                self.status_lbl.config(text='Cannot parse file', fg='#cf222e')
        except Exception as e:
            biosim_log('warn', '_quick_open_file: {}'.format(e))
            self.status_lbl.config(text='Error: ' + str(e)[:50], fg='#cf222e')

    def _parse_sequence_file(self, filepath):
        """Parse FASTA or GenBank file. Returns (sequence, name, type)."""
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # GenBank
        if 'LOCUS' in content[:200] or filepath.lower().endswith(('.gb', '.gbk', '.genbank')):
            try:
                # Extract LOCUS name
                name = 'Unknown'
                for line in content.split('\n'):
                    if line.startswith('LOCUS'):
                        parts = line.split()
                        if len(parts) > 1:
                            name = parts[1]
                        break
                # Extract ORIGIN sequence
                origin_idx = content.find('ORIGIN')
                if origin_idx > 0:
                    seq_part = content[origin_idx:]
                    end_idx = seq_part.find('//')
                    if end_idx > 0:
                        seq_part = seq_part[:end_idx]
                    seq = ''.join(c.upper() for c in seq_part if c.isalpha())
                    return seq, name, 'dna'
            except Exception as _exc:
                biosim_log('warn', '_parse_sequence_file: {}'.format(_exc))
        
        # FASTA
        if content.startswith('>'):
            lines = content.split('\n')
            name = lines[0][1:].split()[0] if len(lines[0]) > 1 else 'Untitled'
            seq = ''.join(c.upper() for c in '\n'.join(lines[1:]) if c.isalpha())
            # Guess type by composition
            if seq and len(set(seq) - set('ACGTNacgtn')) == 0:
                stype = 'dna'
            elif seq and 'U' in seq and 'T' not in seq:
                stype = 'rna'
            else:
                stype = 'protein'
            return seq, name, stype
        
        # Plain text — assume sequence
        seq = ''.join(c.upper() for c in content if c.isalpha())
        name = os.path.basename(filepath).rsplit('.', 1)[0]
        return seq, name, 'dna'

    def _quick_fetch_ncbi(self):
        """Quick NCBI fetch — switch to Databases tab."""
        try:
            self.notebook.select(self.tab_databases)
            self.status_lbl.config(
                text=('Введи accession и жми Fetch' if self.lang == 'ru'
                      else 'Enter accession and click Fetch'),
                fg='#0969da')
        except Exception as _exc:
            biosim_log('warn', '_quick_fetch_ncbi: {}'.format(_exc))

    def _update_status(self):
        """Update sequence info in status bar."""
        if not self.current_sequence:
            self.seq_info_lbl.config(text='')
            return
        seq = self.current_sequence
        info = '{} · {} bp'.format(self.sequence_name[:20], len(seq))
        if self.sequence_type == 'dna':
            gc = (seq.count('G') + seq.count('C')) / len(seq) * 100 if seq else 0
            info += ' · GC: {:.1f}%'.format(gc)
        self.seq_info_lbl.config(text=info)

    # ════════════════════════════════════════════════
    # NEW TABS: Sanger chromatogram, Variant calling, Circular plasmid map
    # ════════════════════════════════════════════════
    def _open_sanger_dialog(self):
        """Sanger sequencing chromatogram viewer (like SnapGene .ab1 viewer).
        Opens as a dialog from the Sequence tab, using the current sequence."""
        is_ru = (self.lang == 'ru')
        tab = tk.Toplevel(self)
        tab.title('Хроматограмма Сэнгера' if is_ru else 'Sanger Chromatogram')
        tab.geometry('820x520')
        tab.configure(bg=biosim_theme('bg'))
        tab.transient(self)
        tk.Label(tab, text=('Хроматограмма Сэнгера' if is_ru
                            else 'Sanger Chromatogram'),
                  bg=biosim_theme('bg'), fg='#0969da', font=('Arial', 13, 'bold')).pack(pady=8)
        tk.Label(tab, text=('Имитация трассы секвенатора для последовательности'
                            if is_ru else 'Simulated sequencer trace for the sequence'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'), font=('Arial', 9)).pack()
        self.sanger_input = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                      font=('Courier', 10), height=3, relief='flat')
        self.sanger_input.pack(fill='x', padx=16, pady=8)
        # Pre-fill with the current sequence from the editor
        cur = self.get_current_sequence()
        self.sanger_input.insert('1.0', cur[:60] if cur else 'ATCGATCGATCGTAGCTAGCTAGCATCGATCG')
        self.sanger_canvas = tk.Canvas(tab, bg=biosim_theme('panel'), height=200,
                                       relief='flat', highlightthickness=1,
                                       highlightbackground=biosim_theme('border'))
        self.sanger_canvas.pack(fill='both', expand=True, padx=16, pady=(0, 6))
        self.sanger_info = tk.Label(tab, text='', bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                     font=('Arial', 9))
        self.sanger_info.pack()
        tk.Button(tab, text=('▶ Показать трассу' if is_ru else '▶ Show trace'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=self._run_sanger).pack(fill='x', padx=16, pady=8)
        tab.after(200, self._run_sanger)

    def _run_sanger(self):
        is_ru = (self.lang == 'ru')
        seq = self.sanger_input.get('1.0', 'end').strip()
        r = simulate_sanger_chromatogram(seq)
        cv = self.sanger_canvas
        cv.delete('all')
        if r.get('error'):
            self.sanger_info.config(text=r['error'], fg='#cf222e')
            return
        cv.update_idletasks()
        W = cv.winfo_width() or 700
        H = 200
        margin = 20
        n = len(r['sequence'])
        if n == 0:
            return
        step = (W - 2 * margin) / max(1, n)
        colors = {'A': '#1a7f37', 'C': '#0969da', 'G': '#1c1c1c', 'T': '#cf222e'}
        baseline = H - 40
        peak_h = baseline - 20
        # Draw peaks for each channel
        for ch in 'ACGT':
            heights = r['channels'][ch]
            pts = []
            for i, h in enumerate(heights):
                x = margin + i * step + step / 2
                y = baseline - h * peak_h
                pts.append((x, y))
            # Draw smooth-ish peaks
            for i in range(len(pts) - 1):
                cv.create_line(pts[i][0], pts[i][1], pts[i+1][0], pts[i+1][1],
                               fill=colors[ch], width=1)
        # Draw base letters
        for i, base in enumerate(r['sequence']):
            x = margin + i * step + step / 2
            if step > 7:
                cv.create_text(x, baseline + 12, text=base, fill=colors.get(base, '#000'),
                               font=('Courier', 8, 'bold'))
        # Mark trim regions (low quality)
        if r['trim_start'] > 0:
            x = margin + r['trim_start'] * step
            cv.create_line(x, 10, x, baseline, fill='#cf222e', dash=(3, 3))
        self.sanger_info.config(
            text=('Средн. качество: {} | Тримминг: {}-{}'.format(
                r['mean_quality'], r['trim_start'], r['trim_end']) if is_ru
                else 'Mean quality: {} | Trim: {}-{}'.format(
                r['mean_quality'], r['trim_start'], r['trim_end'])),
            fg='#1a7f37')

    def _open_variants_dialog(self):
        """Variant calling — compare sample to reference, find SNPs/indels.
        Opens from the Alignment tab, pre-filled with its sequences."""
        is_ru = (self.lang == 'ru')
        tab = tk.Toplevel(self)
        tab.title('Поиск вариантов' if is_ru else 'Variant Calling')
        tab.geometry('760x600')
        tab.configure(bg=biosim_theme('bg'))
        tab.transient(self)
        tk.Label(tab, text=('Поиск вариантов (SNP / indel)' if is_ru
                            else 'Variant Calling (SNP / indel)'),
                  bg=biosim_theme('bg'), fg='#0969da', font=('Arial', 13, 'bold')).pack(pady=8)
        tk.Label(tab, text=('Референс:' if is_ru else 'Reference:'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'), font=('Arial', 9)).pack(anchor='w', padx=16)
        self.var_ref = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                  font=('Courier', 10), height=3, relief='flat')
        self.var_ref.pack(fill='x', padx=16, pady=4)
        # Pre-fill from the alignment tab's sequences if available
        try:
            ref_seq = self.align_seq1_text.get('1.0', 'end').strip()
        except Exception as _exc:
            biosim_log('warn', '_open_variants_dialog: {}'.format(_exc))
            ref_seq = ''
        self.var_ref.insert('1.0', ref_seq or 'ATCGATCGATCGATCGATCGATCG')
        tk.Label(tab, text=('Образец:' if is_ru else 'Sample:'),
                  bg=biosim_theme('bg'), fg=biosim_theme('text_dim'), font=('Arial', 9)).pack(anchor='w', padx=16)
        self.var_sample = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                    font=('Courier', 10), height=3, relief='flat')
        self.var_sample.pack(fill='x', padx=16, pady=4)
        try:
            smp_seq = self.align_seq2_text.get('1.0', 'end').strip()
        except Exception as _exc:
            biosim_log('warn', '_open_variants_dialog: {}'.format(_exc))
            smp_seq = ''
        self.var_sample.insert('1.0', smp_seq or 'ATCGATGGATCGATCGATTGATCG')
        self.var_output = scrolledtext.ScrolledText(tab, bg=biosim_theme('bg'), fg=biosim_theme('text'),
                                                    font=('Courier', 10), height=10, relief='flat', state='disabled')
        self.var_output.pack(fill='both', expand=True, padx=16, pady=8)
        tk.Button(tab, text=('▶ Найти варианты' if is_ru else '▶ Call variants'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=self._run_variants).pack(fill='x', padx=16, pady=8)
        tab.after(200, self._run_variants)

    def _run_variants(self):
        is_ru = (self.lang == 'ru')
        ref = self.var_ref.get('1.0', 'end').strip()
        smp = self.var_sample.get('1.0', 'end').strip()
        r = call_variants(ref, smp)
        out = self.var_output
        out.config(state='normal'); out.delete('1.0', 'end')
        out.tag_config('hdr', foreground='#0969da', font=('Courier', 11, 'bold'))
        out.tag_config('ok', foreground='#1a7f37', font=('Courier', 10, 'bold'))
        out.tag_config('snp', foreground='#cf222e', font=('Courier', 10))
        out.tag_config('indel', foreground='#9a6700', font=('Courier', 10))
        if r.get('error'):
            out.insert('end', r['error'] + '\n'); out.config(state='disabled'); return
        out.insert('end', ('✓ Анализ завершён\n\n' if is_ru else '✓ Analysis complete\n\n'), 'ok')
        out.insert('end', ('Идентичность: {}%\n'.format(r['identity']) if is_ru
                           else 'Identity: {}%\n'.format(r['identity'])), 'hdr')
        out.insert('end', ('Вариантов: {} (SNP: {})\n'.format(r['n_variants'], r['n_snp']) if is_ru
                           else 'Variants: {} (SNP: {})\n'.format(r['n_variants'], r['n_snp'])), 'hdr')
        out.insert('end', ('Ti/Tv: {}\n\n'.format(r['ti_tv_ratio'])), 'hdr')
        for v in r['variants']:
            tag = 'snp' if v['type'] == 'SNP' else 'indel'
            out.insert('end', '  поз {:>5}: {:<12} [{}]\n'.format(
                v['position'], v['change'], v['type']), tag)
        out.config(state='disabled')

    def _open_plasmid_map_dialog(self):
        """Circular plasmid map (like SnapGene's circular view). Opens from the
        Restriction tab, using the current sequence."""
        is_ru = (self.lang == 'ru')
        tab = tk.Toplevel(self)
        tab.title('Карта плазмиды' if is_ru else 'Plasmid Map')
        tab.geometry('700x560')
        tab.configure(bg=biosim_theme('bg'))
        tab.transient(self)
        tk.Label(tab, text=('Кольцевая карта плазмиды' if is_ru
                            else 'Circular Plasmid Map'),
                  bg=biosim_theme('bg'), fg='#0969da', font=('Arial', 13, 'bold')).pack(pady=8)
        self.plasmid_input = scrolledtext.ScrolledText(tab, bg=biosim_theme('panel2'), fg=biosim_theme('text'),
                                                       font=('Courier', 9), height=3, relief='flat')
        self.plasmid_input.pack(fill='x', padx=16, pady=8)
        # Pre-fill with the current sequence
        cur = self.get_current_sequence()
        self.plasmid_input.insert('1.0', cur if cur else
            'ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGATAA' * 3)
        self.plasmid_canvas = tk.Canvas(tab, bg=biosim_theme('panel'), height=320,
                                        relief='flat', highlightthickness=1,
                                        highlightbackground=biosim_theme('border'))
        self.plasmid_canvas.pack(fill='both', expand=True, padx=16, pady=(0, 6))
        tk.Button(tab, text=('▶ Построить карту' if is_ru else '▶ Build map'),
                   bg='#1a7f37', highlightbackground='#1a7f37', fg='white',
                   relief='flat', font=('Arial', 10, 'bold'), cursor='hand2',
                   pady=4, command=self._run_plasmid_map).pack(fill='x', padx=16, pady=8)
        tab.after(200, self._run_plasmid_map)

    def _run_plasmid_map(self):
        import math
        is_ru = (self.lang == 'ru')
        seq = self.plasmid_input.get('1.0', 'end').strip()
        r = build_circular_plasmid_map(seq)
        cv = self.plasmid_canvas
        cv.delete('all')
        if r.get('error'):
            return
        cv.update_idletasks()
        W = cv.winfo_width() or 600
        H = 320
        cx, cy = W / 2, H / 2
        radius = min(W, H) / 2 - 50
        # Draw the main circle (plasmid backbone)
        cv.create_oval(cx - radius, cy - radius, cx + radius, cy + radius,
                       outline='#1c2128', width=3)
        # Center label
        cv.create_text(cx, cy - 8, text=('Плазмида' if is_ru else 'Plasmid'),
                       fill='#1c2128', font=('Arial', 11, 'bold'))
        cv.create_text(cx, cy + 10, text='{} bp'.format(r['length']),
                       fill='#57606a', font=('Arial', 10))
        cv.create_text(cx, cy + 26, text='GC {}%'.format(r['gc_percent']),
                       fill='#9a6700', font=('Arial', 9))
        # Draw features as colored arcs
        feat_colors = ['#1a7f37', '#0969da', '#cf222e', '#9a6700', '#8250df', '#bc4c00']
        for i, f in enumerate(r['features']):
            a1 = f['start_angle']; a2 = f['end_angle']
            if a2 <= a1:
                a2 = a1 + 5
            col = feat_colors[i % len(feat_colors)]
            # Arc on the circle (tkinter: 0°=east, counterclockwise; convert)
            start = 90 - a2
            extent = a2 - a1
            cv.create_arc(cx - radius, cy - radius, cx + radius, cy + radius,
                          start=start, extent=extent, style='arc',
                          outline=col, width=8)
            # Feature label at midpoint angle
            mid = math.radians(90 - (a1 + a2) / 2)
            lx = cx + (radius + 22) * math.cos(mid)
            ly = cy - (radius + 22) * math.sin(mid)
            cv.create_text(lx, ly, text=f['name'], fill=col, font=('Arial', 8, 'bold'))
        # Position 0 marker (top)
        cv.create_line(cx, cy - radius - 8, cx, cy - radius + 8, fill='#cf222e', width=2)
        cv.create_text(cx + 14, cy - radius - 4, text='0', fill='#cf222e', font=('Arial', 8))



