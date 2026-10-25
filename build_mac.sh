#!/bin/bash
# ================================================================
#  VIREON — сборка VIREON.app для macOS
#  Запуск в Терминале:   bash build_mac.sh
#  Готовый .app появится в папке  dist/
# ================================================================
set -e
cd "$(dirname "$0")"

echo
echo "=== [1/5] Проверяю Python ==="
python3 --version || { echo "ОШИБКА: нет python3. Поставь с https://www.python.org/downloads/"; exit 1; }

echo
echo "=== [2/5] Ставлю обязательные библиотеки ==="
python3 -m pip install --upgrade pip
python3 -m pip install numpy scipy matplotlib Pillow
python3 -m pip install pyinstaller

echo
echo "=== [3/5] Чищу прошлую сборку и кэш PyInstaller ==="
rm -rf build dist VIREON.spec
rm -rf "$HOME/Library/Application Support/pyinstaller"

# Обход ошибки 'Executable contains code signature!' на старых macOS
# (баг самого PyInstaller, его issue #6167). Патчим установленный
# PyInstaller так, чтобы он снимал лишнюю подпись, а не падал.
echo "--- применяю обход бага подписи PyInstaller ---"
python3 fix_pyinstaller.py || echo "(патч уже применён или не требуется)"

echo
echo "=== [4/5] Собираю VIREON.app (пара минут) ==="
# Тяжёлые опциональные библиотеки ИСКЛЮЧАЕМ из сборки:
# программе они для запуска не нужны (грузятся мягко через try/except),
# а rdkit/torch и прочие ломают сборку на старом macOS.
python3 -m PyInstaller --windowed --name VIREON \
    --add-data "body_organs.png:." \
    --paths "." \
    --hidden-import biosim_part2 \
    --hidden-import biosim_part3 \
    --hidden-import biosim_part4 \
    --collect-submodules matplotlib \
    --exclude-module rdkit \
    --exclude-module torch \
    --exclude-module tensorflow \
    --exclude-module sklearn \
    --exclude-module scikit-learn \
    --exclude-module deepchem \
    --exclude-module esm \
    --exclude-module cobra \
    --exclude-module libsbml \
    --exclude-module sbol2 \
    --exclude-module py3Dmol \
    --exclude-module pydicom \
    --exclude-module h5py \
    --exclude-module plotly \
    --exclude-module statsmodels \
    --exclude-module nbformat \
    --exclude-module pandas \
    --exclude-module Bio \
    --exclude-module sympy \
    --exclude-module PyQt5 \
    --exclude-module PySide2 \
    biology.py

echo
echo "=== [5/5] Готово ==="
echo "Твоя программа:  dist/VIREON.app"
echo "Запусти так:  open dist/VIREON.app"
echo "Можешь передавать её людям — Python им НЕ нужен."
