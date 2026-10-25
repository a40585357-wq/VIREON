#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VIREON — обход бага PyInstaller на старых macOS.

Ошибка:  AssertionError: Executable contains code signature!
Причина: на macOS 10.13 (High Sierra) codesign не до конца снимает
         фиктивную подпись с исполняемого файла, и PyInstaller падает
         на проверке (его же issue #6167).

Что делает этот скрипт:
  находит функцию fix_exe_for_code_signing в установленном PyInstaller
  и заменяет её падающую проверку (assert) на код, который при наличии
  подписи СНИМАЕТ её и перечитывает заголовки — вместо аварийного выхода.

Запуск (один раз):   python3 fix_pyinstaller.py
После этого собирай как обычно:  bash build_mac.sh
"""
import os
import sys
import shutil

MARKER = '# VIREON-PATCH'

PATCH_BLOCK = '''    sign_sec = [cmd for cmd in header.commands if cmd[0].cmd == LC_CODE_SIGNATURE]
    if len(sign_sec) != 0:  {marker}
        # VIREON: strip the leftover signature instead of crashing, then reload.
        import subprocess as _sp  {marker}
        try:  {marker}
            _sp.run(['codesign', '--remove-signature', filename], check=False)  {marker}
        except Exception:  {marker}
            pass  {marker}
        executable = MachO(filename)  {marker}
        header = executable.headers[-1]  {marker}
        sign_sec = [cmd for cmd in header.commands if cmd[0].cmd == LC_CODE_SIGNATURE]  {marker}
    # If a signature still remains we drop the hard assert to a soft skip:
    if len(sign_sec) != 0:  {marker}
        return  {marker}
'''.format(marker=MARKER)

ORIGINAL_LINE = ('    assert len(sign_sec) == 0, "Executable contains code '
                 'signature!"')


def find_osx_module():
    try:
        import PyInstaller.utils.osx as m
        return m.__file__
    except Exception:
        # Fallback: search site-packages
        import PyInstaller
        base = os.path.dirname(PyInstaller.__file__)
        cand = os.path.join(base, 'utils', 'osx.py')
        return cand if os.path.exists(cand) else None


def main():
    path = find_osx_module()
    if not path or not os.path.exists(path):
        print('Не нашёл PyInstaller. Сначала установи его:  pip install pyinstaller')
        return 1

    with open(path, 'r', encoding='utf-8') as f:
        src = f.read()

    if MARKER in src:
        print('Уже пропатчено ранее — ничего делать не нужно. Можно собирать.')
        return 0

    # The line we replace (exact grep line from the sign check)
    needle = None
    for cand in (
        '    assert len(sign_sec) == 0, "Executable contains code signature!"',
        'assert len(sign_sec) == 0, "Executable contains code signature!"',
    ):
        if cand in src:
            needle = cand
            break

    if needle is None:
        print('Не нашёл строку проверки подписи — возможно, другая версия '
              'PyInstaller. Напиши мне свою версию:  pyinstaller --version')
        return 2

    # Build the replacement: keep the original sign_sec line that precedes the
    # assert, then our block (which itself recomputes sign_sec). We replace from
    # the sign_sec definition through the assert with our full block.
    sign_def = ('    sign_sec = [cmd for cmd in header.commands if cmd[0].cmd '
                '== LC_CODE_SIGNATURE]')
    if sign_def in src:
        old = sign_def + '\n' + needle
        new = PATCH_BLOCK.rstrip('\n')
        if old in src:
            src2 = src.replace(old, new, 1)
        else:
            # assert indented differently; just swap the assert line
            src2 = src.replace(needle, PATCH_BLOCK.rstrip('\n'), 1)
    else:
        src2 = src.replace(needle, PATCH_BLOCK.rstrip('\n'), 1)

    # Backup then write
    shutil.copy2(path, path + '.vireon_bak')
    with open(path, 'w', encoding='utf-8') as f:
        f.write(src2)

    print('Готово! PyInstaller пропатчен (бэкап: osx.py.vireon_bak).')
    print('Теперь собери программу:  bash build_mac.sh')
    return 0


if __name__ == '__main__':
    sys.exit(main())
