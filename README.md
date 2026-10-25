# VIREON

*English below · [Русская версия ниже](#vireon--русская-версия)*

**Cellular & Synthetic Biology Suite** — a desktop simulator for cellular and
synthetic biology.

VIREON brings genetic-circuit design, genome editing, metabolic modeling,
whole-cell simulation, and molecular-biology tools together in a single program.

---

## Features

- **Constructor** — a visual designer for synthetic genetic circuits:
  parts on a continuous DNA backbone (SBOL glyphs), reaction templates
  (1→1, 2→1, etc.), regulation (activation / repression / transcription /
  direct gene regulation), a Flux Heatmap, compartmentalization (organelles),
  an ATP energy currency, deterministic ODE runs, canvas zoom and scrolling.
- **CRISPR** — guide-RNA selection and analysis.
- **Metabolism** — pathway models (ODE / RK4, Gillespie, FBA/FVA/pFBA),
  linked to metabolic diseases (e.g. diabetes → disrupted pathway).
- **Synthetic organism** — whole-cell simulation.
- **Molecular biology** — sequence handling, drug analysis and drug design.
- **Main window** — human anatomy, pharmacology, disease selection, and a
  detailed disease picture (genes, proteins, cascade) with sequences pulled
  from NCBI.
- Light **Ocean** theme and dark **Graphite Neon** theme (toggled with the
  day/night button).

---

## Install & run

Requires **Python 3.8+** (the `tkinter` module ships with Python).

```bash
# 1. install dependencies
pip install -r requirements.txt

# 2. run
python biology.py
```

> `biology.py` is the entry point; the other three files (`biosim_part2.py`,
> `biosim_part3.py`, `biosim_part4.py`) load automatically and must sit in the
> same folder.

Optional extensions (COBRA, RDKit, BioPython, PyTorch, etc.) are installed only
if you want them — see the comments in `requirements.txt`. The core works
without them.

---

## Project layout

| File               | Purpose                                                   |
|--------------------|-----------------------------------------------------------|
| `biology.py`       | Main window, entry point, anatomy/pharmacology            |
| `biosim_part2.py`  | Synthetic-biology window (Constructor, CRISPR, etc.)      |
| `biosim_part3.py`  | Computation engines (metabolism, organelles, ATP)         |
| `biosim_part4.py`  | Molecular-biology window                                  |

---

## Acknowledgments & sources of ideas

The logic of several sections is inspired by open scientific tools. Their code
was **not copied** — only ideas and principles were borrowed, rethought, and
reimplemented in Python:

- **TinkerCell** — approach to visual genetic-circuit design.
- **COPASI** — metabolic/kinetic-modeling principles.
- **CRISPOR** — guide-RNA selection logic.

---

## License

See the [`LICENSE`](LICENSE) file (proprietary license).

In short:
- **Free** — to download and **use** the program for personal, educational, and
  other **noncommercial** purposes.
- **Not allowed** — to modify the program or create other programs or
  derivative products based on it.
- **Commercial use** — only under a separate **paid** license (monthly rental).
  Contact the author.

===========================================================================

# VIREON — Русская версия

*[English version above](#vireon)*

**Cellular & Synthetic Biology Suite** — настольный симулятор клеточной и
синтетической биологии.

VIREON объединяет в одной программе проектирование генетических схем,
редактирование генома, метаболическое моделирование, моделирование целой
клетки и инструменты молекулярной биологии.

---

## Возможности

- **Конструктор** — визуальный дизайнер синтетических генетических схем:
  части на непрерывном ДНК-остове (SBOL-глифы), шаблоны реакций (1→1, 2→1 и др.),
  регуляция (активация / репрессия / транскрипция / прямая генная регуляция),
  тепловая карта потоков (Flux Heatmap), компартментализация (органеллы),
  энергетическая валюта ATP, детерминированный ODE-запуск, зум и прокрутка холста.
- **CRISPR** — подбор и анализ направляющих РНК.
- **Метаболизм** — модели путей (ODE / RK4, Gillespie, FBA/FVA/pFBA),
  связь с метаболическими болезнями (например, диабет → нарушенный путь).
- **Синтетический организм** — моделирование целой клетки.
- **Молекулярная биология** — работа с последовательностями, анализ и создание лекарств.
- **Главное окно** — анатомия человека, фармакология, выбор болезни и
  построение детальной картины болезни (гены, белки, каскад) с загрузкой
  последовательностей из NCBI.
- Светлая тема **Ocean** и тёмная тема **Graphite Neon** (переключение кнопкой день/ночь).

---

## Установка и запуск

Нужен **Python 3.8+** (модуль `tkinter` входит в стандартную поставку Python).

```bash
# 1. установить зависимости
pip install -r requirements.txt

# 2. запустить
python biology.py
```

> `biology.py` — точка входа; три остальных файла (`biosim_part2.py`,
> `biosim_part3.py`, `biosim_part4.py`) загружаются автоматически и должны
> лежать в той же папке.

Необязательные расширения (COBRA, RDKit, BioPython, PyTorch и др.) ставятся
по желанию — см. комментарии в `requirements.txt`. Без них работает ядро программы.

---

## Структура проекта

| Файл              | Назначение                                             |
|-------------------|--------------------------------------------------------|
| `biology.py`      | Главное окно, точка входа, анатомия/фармакология        |
| `biosim_part2.py` | Окно синтетической биологии (Конструктор, CRISPR и др.) |
| `biosim_part3.py` | Расчётные движки (метаболизм, органеллы, ATP)           |
| `biosim_part4.py` | Окно молекулярной биологии                              |

---

## Благодарности и источники идей

Логика ряда разделов вдохновлена открытыми научными инструментами. Их код
**не копировался** — заимствовались только идеи и принципы, переосмысленные и
реализованные на Python:

- **TinkerCell** — подход к визуальному дизайну генетических схем.
- **COPASI** — принципы метаболического/кинетического моделирования.
- **CRISPOR** — логика подбора направляющих РНК.

---

## Лицензия

См. файл [`LICENSE`](LICENSE) (проприетарная лицензия).

Коротко:
- **Бесплатно** — скачивать и **пользоваться** программой в личных, учебных
  и других **некоммерческих** целях.
- **Запрещено** — изменять программу и создавать на её основе другие
  программы или производные продукты.
- **Коммерческое использование** — только по отдельной **платной** лицензии
  (ежемесячная аренда). По вопросам — свяжитесь с автором.
