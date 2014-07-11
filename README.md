# *catm* - cellular automata

This repository is a collection of functions and code usable to run cellular automata simulations in R. The aim is to develop it into a valid and distributable R-package. 
 
For now, it contains standalone files with examples that were published in peer-review articles:

  - Conway's game of life: [version1](/gameoflife.r). This is the most famous cellular automaton. Read more about it on [Wikipedia.org](https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life). 
  - A vegetation facilitation model (Kéfi et al. 2008): [version1](/facilitationmodel.r).  This is Sonias application of cellular automata to simulate spatial vegetation structure from simple assumptions on facilitation. The cellular automata used in cascade are building upon this model. 



## R

You will need R to run the code. Download the latest version for your operating system from [here](http://stat.ethz.ch/CRAN/) and install it. It is a command line based program. You just copy phrases or lines of code into the R terminal and they immediately return a value. For more convenient use of R, I recommend using an integrated development environment (IDE) such as RStudio. 

## RStudio: an IDE for R

You could manipulate the R code  with any text editor and copy the lines of code to the R terminal. Chose whatever you prefer to open and manipulate the files. However, for R there is a very convenient programming environment called [RStudio](https://www.rstudio.com/ide/download/desktop). It provides you with all you need to use R, like the documentation of functions, a file browser and a preview of the available R objects. 

## Animated output with *ImageMagick*

If you like to create animated gifs as a visualisation of the simulation you will need to install another little program called ImageMagick. You can download it [here](http://www.imagemagick.org/script/binary-releases.php).  Install it whereever you like, but you need to make sure


## License

The MIT License (MIT)

Copyright (c) 2014 Florian D. Schneider

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.