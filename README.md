# EvaalScore manual

2024-11-8

This manual is a work in progress.

## Preparing a working directory for EvaalScore

You run `evaalscore` from a working directory.  Create a `project` directory under it and a file called `distance` in the `project` directory.  The file `distance` contains one line of text: either "localdistance" or "wgs84distance", depending on the reference system used for the coordinates.  It must be either a local system with unit in meters or else WGS84. In the latter case, the coordinates must be in the order lon,lat.

You need a `groundTruth.txt` file to get a minimum sets of statistics.  Using this file, EvaalScore displays some path characteristics, like total lenght and number of points per floor.

The files `buttonsPressed.log` and `positions.log` are created by the stepLogger app.  When those files are present, EvaalScore makes its whole sets of statistics and prints them.  It additionaly creates output files `qgis_errors.txt` and `error.report`.

### Creating a dummy path for testing


## Working with image maps

If you plan to use maps to obtain a graphic report, you need to put a file called `legend.png` inside the `project` directory.  You also need to have georeferenced maps and to create a scenario describing the positioning of the maps, legend and statistics on the graphic report.


### Creating `image.png` and `image.png.points` with the Qgis georeferenceer

Maps are used as background images on the output of `evaalsore`, that is files `path.pdf`, `score.pdf` and `path.ogg`.  All of these need maps as background on which both the reference points and the estimated path are drawn.

You need PNG maps with well-visible lines.  In practice, for Qgis this means that the PNG file size should generally be around 500kB with around 4000 pixels on the widest dimension.  Smaller images risk having too ragged lines, while larger will have too thin lines.  For EvaalScore you need smaller images, generally around 100kB in size with around 1000 pixel along the widest dimension.

You should start from a high quality source, ideally a vector map, and convert to a high-density PNG image, for example with 12k pixels on the widest dimension.  Then run the `shrink-png` bash script to reduce the size of this initial image to what is appropriate for Qgis.  Experiment with the initial PNG image size and the number of iterations needed to reduce the image size to the ideal size.  `shrink-png` uses a default number of 4 iterations, which means shrinking to 60%: look at the comment in the source for how to change the number of iterazions.

Once you have loaded the image into the Qgis georeferencer and having georeferenced it, you can produce an `image.png.points` file. Copy it together with `image.png` into the `project` directory.  Repeat the procedure for all maps you need for this project.  Then, from the working directory, run `evaalscore` and choose Exit from the menu: this defines the `create_corners` function.  Now, for each `image.png` file, run `create_corners image project`: this creates the `image.png.corners` file.  If for some reason you have no `image.png.points` file, you can create the `image.png.corners` file by hand: look at the comments inside `evaalscore.m` for the format of the file.


### Creating a scenario

### Preparing the corners file and shrinking the final map images

Use create_corners to create an `image.corners` for each `image.png`, which creates an `image.corners` file with the lon,lat coordinate of the four corners of the image.  After this is done, the `image.png`points` file is no more needed, and you can shrink the images to make thin lines more visible.

Shring all the maps using `shrink-png` so that lines are more clearly visible in the small images that EvaalScore uses for its output: `path.pdf`, `score.pdf` and `path.ogg` all have maps as background on which both the reference points and the estimated path are drawn.


### Creating dummy data files for testing the graphical output of a scenario

If you have the map images and you have defined a scenario you miss the data files, you can create dummy ones for testing how the graphics output will look like.

You can create a dummy `groundTruth.txt` file using the following code, after having read the first two columns from image.corners into the corners variable

	pkg load matgeom
	corners = [
	  114.178206,22.302603
	  114.178678,22.304240
	  114.179900,22.303930
	  114.179428,22.302293];
    minc=min(corners);
	ranc=range(corners);
	floor=3;
	N = 30;
	refp = zeros(N,4);
	do
	  refp(1,:) = [1, rand(1,2).*ranc+minc, floor];
	until(isPointInPolygon(refp(1,2:3), corners))
	for p = 2:N
	  do
	    refp(p,:) = [p, refp(p-1,2:3)+randn(1,2).*ranc*10/N, floor];
	  until(isPointInPolygon(refp(p,2:3), corners))
	endfor
	fid = fopen("groundTruth.txt", 'w');
	fprintf(fid, "%3d%12.6f%12.6f%3d\n", refp');
	fclose(fid);

If you have a `groundTruth.txt` file but either `buttonsPressed.log` or `positions.log` do not exist in the project directory, the Evalscore menu does not offer to use that directory.  However, if `dirname` is set and `dopath==true`, when calling `evaalscore`, it uses `groundTruth.txt` as a starting point to create dummy `buttonsPressed.log` and `positions.log` files in their place so that a `path.pdf` file can be created.  In this case, statistics that EvaalScore prints out and includes in the `score.pdf` report make no sense as the dummy files created emulate an error-free run.

It is also possible to generate a dummy path with errors for testing how `score.pdf` will look like.  Calling `create_dummy_positions(dirnanme, 5)` will create `buttonsPressed.log` and `positions.log` files with random errors and a `score.pdf` report with meaningful yet dummy statistics.
