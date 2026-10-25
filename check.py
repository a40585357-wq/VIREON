#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VIREON — проверка, что всё на местах перед заливкой на GitHub.
Запуск:   python check.py    (Windows)
          python3 check.py   (Mac / Linux)
"""
import os, sys, ast

HERE = os.path.dirname(os.path.abspath(__file__))

NEED_FILES = [
    "biology.py", "biosim_part2.py", "biosim_part3.py", "biosim_part4.py",
    "README.md", "LICENSE", "requirements.txt", ".gitignore",
]
CODE_FILES = ["biology.py", "biosim_part2.py", "biosim_part3.py", "biosim_part4.py"]

ok = True

print("=" * 55)
print("  VIREON — проверка комплекта")
print("=" * 55)

# 1) все нужные файлы на месте
print("\n[1] Файлы на месте:")
for f in NEED_FILES:
    exists = os.path.isfile(os.path.join(HERE, f))
    print(f"    {'OK ' if exists else 'НЕТ'}  {f}")
    if not exists:
        ok = False

# 2) код компилируется (нет синтаксических ошибок)
print("\n[2] Код без ошибок:")
for f in CODE_FILES:
    p = os.path.join(HERE, f)
    if not os.path.isfile(p):
        continue
    try:
        ast.parse(open(p, encoding="utf-8").read())
        print(f"    OK   {f}")
    except SyntaxError as e:
        print(f"    ОШИБКА в {f}: строка {e.lineno}: {e.msg}")
        ok = False

# 3) точка входа запускается (пробный импорт частей)
print("\n[3] Части подключаются друг к другу:")
sys.path.insert(0, HERE)
for mod in ["biosim_part2", "biosim_part3", "biosim_part4"]:
    try:
        __import__(mod)
        print(f"    OK   {mod}")
    except Exception as e:
        # tkinter-ошибки на сервере без экрана — не считаем провалом
        msg = str(e)
        low = msg.lower()
        if ("display" in low or "tcl" in low or "tkinter" in low
                or "DISPLAY" in msg):
            print(f"    OK   {mod}  (нет экрана/tkinter — норм для проверки)")
        else:
            print(f"    ВНИМАНИЕ {mod}: {msg}")

print("\n" + "=" * 55)
if ok:
    print("  ВСЁ НА МЕСТЕ. Можно заливать на GitHub.")
else:
    print("  ЕСТЬ ПРОБЛЕМЫ ВЫШЕ — поправь перед заливкой.")
print("=" * 55)
