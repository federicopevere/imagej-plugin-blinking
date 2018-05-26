# imagej-plugin-blinking

This is my repository in Git about an ImageJ plugin which can be used to extract
intensity time traces from intensity image sequences of luminescent
semiconductor quantum dots.
Brief description of how the plugin works:
1. After installing the plugin, open ImageJ and load the image sequence in TIFF
   format: File -> Import -> Image sequence . To avoid very slow processing of
   the image sequence it is important not to load it as "Virtual stack".
2. The image sequence must be stabilized before use. In order to do that, one
    use for example Kang Li’s Image Stabilizer.
3. Click on Images - > Stacks -> Z project: select Average and save the image
   which will be used to find the dot xy coordinates (also called dot map).
4. Open the image obtained from 3., select your region of interest (ROI) by
   using polygon selection: this area is where your dots are located.
5. Click on Process -> Find Maxima : this will allow to find the position of the
   dots. Select "Preview point selection" adjust the "Noise tolerance" in such a
   way the your dots are marked. Select “List” as output type and save it.
6. For the blinking plugin to work, the object map must have at each row the x
   and y coordinates of one dot separated by a tab character. Thus you need to
   modify the dot map obtained at 5. by using Excel or similar software.
7. Close the image obtained at 3. and now you will be left with only the image
   sequence opened.
8. Run the blinking plugin
   - Select the dot map file
   - Choose not to input a background map
   - Choose the dot edge (must be an odd number) and the background distance
     from the dot edge
   - Input the frame time in seconds
   - Save the output files: total intensity, background and signal time traces
   - Select what traces you want to plot for all the dots


Copyright: Federico Pevere and Benjamin Bruhn (2018), pevere@kth.se

Please put the following references in your publication if you use this plugin for your work:

Text References:

Pevere, F. , Bruhn, B. , Sangghaleh, F. , Hormozan, Y. , Sychugov, I. and Linnros, J. (2015), Effect of X‐ray irradiation on the blinking of single silicon nanocrystals. Phys. Status Solidi A, 212: 2692-2695.

Bibtex:

@article{blinkingSiQDs,
 author = {Pevere Federico and Bruhn Benjamin and Sangghaleh Fatemeh and Hormozan Yashar and Sychugov Ilya and Linnros Jan},
 title = {Effect of X‐ray irradiation on the blinking of single silicon nanocrystals},
 journal = {physica status solidi (a)},
 volume = {212},
 number = {12},
 pages = {2692-2695},
 doi = {10.1002/pssa.201532652}
 }
