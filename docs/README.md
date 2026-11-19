<!--
SPDX-FileCopyrightText: 2015 UCLA

SPDX-License-Identifier: BSD-3-Clause
-->

# VERB4D Doxygen documentation
To generate or update documentation, run the following:
`doxygen Doxyfile`

## Sphinx documentation (pydata-sphinx-theme)
The same theme as verb4dpy is available through Breathe and Exhale. Doxygen must be installed.

```
uv run --no-project --with-requirements requirements.txt sphinx-build -b html . _build/html
```
