c-ray-fast
==========

An improved, but still extremely simple, raytracer.

This version of c-ray is cherry-picking the best parts of c-ray-f and c-ray-mt
which have diverged over time, and adds some sorely needed improvements, while
still maintaining the overall simplicity of c-ray and the simple 0-dependency
build.

Main improvements:
  - simple ray intersection acceleration
  - better task parallelization than c-ray-mt
  - improved backwards-compatible scene format
  - option for specifying the background color
  - option for producing an alpha mask
  - colored lights
  - gamma correction

License
-------
Copyright (C) 2006-2021 John Tsiombikas <nuclear@member.fsf.org>

c-ray-fast is free software. Feel free to use, modify, and/or redistribute it
under the terms of the GNU General Public License v3, or at your option any
later version published by the Free Software Foundation. See COPYING for
details.
