@echo off
REM ================================================================
REM  VIREON - сборка VIREON.exe для Windows
REM  Двойной клик по этому файлу ИЛИ запуск из терминала.
REM  Готовый .exe появится в папке  dist\
REM ================================================================

echo.
echo === [1/5] Проверяю Python ===
python --version
if errorlevel 1 (
    echo.
    echo ОШИБКА: Python не установлен или не в PATH.
    echo Скачай с https://www.python.org/downloads/  и при установке
    echo обязательно поставь галочку "Add Python to PATH".
    pause
    exit /b 1
)

echo.
echo === [2/5] Ставлю обязательные библиотеки ===
python -m pip install --upgrade pip
python -m pip install numpy scipy matplotlib Pillow
python -m pip install pyinstaller

echo.
echo === [3/5] Чищу прошлую сборку ===
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist VIREON.spec del /q VIREON.spec

echo.
echo === [4/5] Собираю VIREON.exe (это займёт пару минут) ===
REM Тяжёлые опциональные библиотеки исключаем из сборки.
pyinstaller --onefile --windowed --name VIREON ^
    --add-data "body_organs.png;." ^
    --paths "." ^
    --hidden-import biosim_part2 ^
    --hidden-import biosim_part3 ^
    --hidden-import biosim_part4 ^
    --collect-submodules matplotlib ^
    --exclude-module rdkit ^
    --exclude-module torch ^
    --exclude-module tensorflow ^
    --exclude-module sklearn ^
    --exclude-module scikit-learn ^
    --exclude-module deepchem ^
    --exclude-module esm ^
    --exclude-module cobra ^
    --exclude-module libsbml ^
    --exclude-module sbol2 ^
    --exclude-module py3Dmol ^
    --exclude-module pydicom ^
    --exclude-module h5py ^
    --exclude-module plotly ^
    --exclude-module statsmodels ^
    --exclude-module nbformat ^
    --exclude-module pandas ^
    --exclude-module Bio ^
    --exclude-module sympy ^
    --exclude-module PyQt5 ^
    --exclude-module PySide2 ^
    biology.py

echo.
echo === [5/5] Готово ===
echo Твой файл лежит здесь:  dist\VIREON.exe
echo Можешь передавать его людям - Python им НЕ нужен.
echo.
pause
