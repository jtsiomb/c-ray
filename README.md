c-ray archival repo
===================

![shots](http://nuclear.mutantstargoat.com/sw/c-ray/img/c-ray_thumbs.jpg)

Quick directions
----------------
If you're looking for the version of c-ray widely used as a benchmark, you can
find it under the `c-ray-mt` directory, by grabbing the revision tagged as
"c-ray-1.1" (it has accrued random hacks in later revisions). You can also get
the original release tarball, exactly as it was packaged and distributed by Ian
Mapleson from any of these mirrors:

  - https://github.com/jtsiomb/c-ray/releases/download/c-ray-1.1/c-ray-1.1.tar.gz
  - http://nuclear.mutantstargoat.com/sw/c-ray/c-ray-1.1.tar.gz

If you're looking for the improved version of the c-ray code (bugs fixed, much
faster, gamma corrected output, more flexible scene description, but still
maintaining the simplicity and zero dependencies of the original), look in the
`c-ray-fast` directory. Or grab the last release archive:

  - http://github.com/jtsiomb/c-ray/releases/download/c-ray-2.0/c-ray-2.0.tar.gz
  - http://nuclear.mutantstargoat.com/sw/c-ray/c-ray-2.0.tar.gz

As a renderer c-ray-2.0 is a *significant* improvement over the original. But it
may or may not be better as a CPU benchmark. Don't rush to upgrade your
benchmarks blindly.

License
-------
Copyright (C) 2005-2021 John Tsiombikas <nuclear@member.fsf.org>

Everything in this repository is free software. Feel free to use, modify, and/or
redistribute under the terms of the GNU General Public License v3, or at your
option any later version published by the Free Software Foundation. See COPYING
for details.

History
-------
I wrote the first version of c-ray one afternoon in 2005, because I was tired of
the much more elaborate C++ renderer I was writing at the time, and wanted to
see how many lines of C code it would take to write an extremely simple
self-contained raytracer in a single source file. It rendered spheres, had
reflections and a line-based scene description, and used SDL 1.2 to continuously
draw frames in a window, featuring a rotating camera and a bouncing ball.

## c-ray-f
I quickly realized that using SDL to draw in a window runs counter to the idea
of making the renderer as simple and short as possible. So I removed the SDL
code, and made it into a filter (the `-f` suffix stands for filter) which reads a
scene description from stdin, and writes an image to stdout, minimizing the code
further and removing any library dependencies other than libc in the process.

There was also an intermediate `c-ray-unix` version, which was similar to the
original SDL one, but using Xlib to talk to the window system instead of SDL.

## Benchmarking PCs, the Gameboy Advance, and SGI workstations
Since `c-ray-f` was so simple to compile and run on any system with a C
compiler, and since it prints the elapsed time to stderr, it turned out to be
very convenient for comparing pure floating point performance.

I ran it on all my computers, including my Gameboy
Advance (you can find the GBA port under `c-ray-gba`), and my SGI Indy, and I've
also sent it to a number of friends, collecting the results in a web page:
http://nuclear.mutantstargoat.com/sw/c-ray/c-ray_results.html

At about this time I was talking to Ian Mapleson, legendary 2nd hand SGI
refurbisher/reseller, because I intended to buy an Octane2 workstation. I
noticed on his website that he had a benchmarking section comparing the
performance of various SGI workstations, but was lacking anything that can test
floating point performance. So I pointed c-ray to him:

> P.S.2 I noticed you are interested in benchmarking, and that you don't          
> have any floating point benchmarks there. I happen to have a simple ray         
> tracing program I wrote for benchmarking reasons. I even used it on             
> various machines recently including an O2 and my Indy. The results,             
> along with the source code, are located here:                                   
> http://apricot.hep.ntua.gr/~nuclear/tmp/c-ray_results.html                      
> if you find it useful, feel free to use it.

## Multithreaded c-ray-mt
Back in 2005 multi-core personal computers were still quite rare, but since Ian
used to deal with high-end SGI systems which commonly came in multi-processor
configurations, it quickly became apparent that the single-threaded `c-ray-f` is
of limited use when it comes to benchmarking multi-CPU or multi-core systems.
Therefore I quickly hacked a simple multi-threaded version of c-ray, called
`c-ray-mt`. The last versions of `c-ray-mt` in the repo accrued a number of
oddball additions and experiments over the years; the original version is the
one tagged with `c-ray-1.1` (or you could just get the release archive).

## Spread within the SGI community
Ian took over the benchmark results page, and it quickly expanded to become a
comprehensive comparison of a slew of different SGI systems:
https://web.archive.org/web/20200209145437/http://www.futuretech.blinkenlights.nl/c-ray.html

Ian's c-ray results page popularized c-ray to the SGI enthusiast community, back
then mostly centered around the now-defunct nekochan forums. The denizens of
nekochan were all too eager to show off their SGI system's performance, and as a
result c-ray spread through the community.

## Spreading far and wide
Seeing c-ray used by nekochan users was interesting, but not very surprising,
given the Ian connection. But then suddenly I started hearing about c-ray
popping up all over the place. Apparently Michael Larabel of Phoronix (a UNIX
hardware review and news site) discovered c-ray, presumably from nekochan, and
started using it for his CPU reviews. It was added to the Phoronix Test Suite,
and became part of OpenBenchmarking.org, at which point it got picked up and
used by more mainstream hardware review sites, like Anandtech and Tom's
Hardware.

It was quite amusing for me because at the time I was working for
Maxon on Cinema4D, and I kept seeing hardware review sites present results
from both c-ray and CineBench (a version of Cinema4D used as a benchmark) one
after the other on the same page.

Finally a few years ago I even saw the CEO of AMD, on stage for an AMD product
launch, referring to c-ray as an "industry-standard benchmark", and using it to
show how much faster is their new processor than the competition:
https://www.youtube.com/watch?v=TtdFOdj_kCQ

C-ray versions
--------------
  - `c-ray.c`: very first version, animated real-time, using SDL for video output.
  - `c-ray-unix`: similar to `c-ray.c`, but dropping SDL in favour of Xlib.
  - `c-ray-f`: minimal filter version, usable in a pipeline. Outputs a ppm
    image, instead of displaying the result on screen.
  - `c-ray-mt`: multi-threaded version of `c-ray-f`. Later versions became a
    playground for random experiments, getting a Torrance&Cook reflectance
    model, and very rudimentary photon mapping for caustics.
  - `c-ray-gba`: Gameboy Advance port.
  - `c-ray-gpu`: animated real-time version running on the GPU as a GLSL shader.
  - `cray16`: 16-bit MS-DOS port (Turbo C), rendering in 256 color VGA modes
    (mode 13h and mode-x) with optional floyd-steinberg dithering.
  - `cray32`: 32-bit MS-DOS port (Watcom), using VESA VBE for high resolution
    and high color/true color modes.
  - `c-ray-fast`: new improved multi-threaded c-ray, with an octree for ray
    intersection acceleration, gamma-correct output, more flexible scene
    description format, and other fixes and improvements.
