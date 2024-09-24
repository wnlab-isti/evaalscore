## Score computer for the EvAAL competition
##
## Copyright © Francesco Potortì, 2015-2016
## Licensed to anyone under the GNU Affero General Public License 3.0

## $Revision: 3.18 $

## When this script is called from Octave it looks under BASEDIR for
## directories containing all the necessary files (see below) and shows
## them in a menu for the user to choose from.  After the choice is
## done, data are read and a report is created.  If doplot is true
## a graph displaying the estimated path superimposed over the map
## is displayed and saved under the chosen dir in a pdf file.

## The normal way to call this script from the Octave prompt is:
## >> evaalscore
## which creates a score.pdf file or, alternatively,
## >> domovie=true evaalscore
## which creates a movie.ogg file instead.

## Setting
## >> dopath=true
## additionally creates a path.pdf file.

## The menu shows the usable subdirectories of basedir, defaulting to "."
## >> basedir="../competitors/data"; evaalscore
## If you want to completely avoid the prompt asking to choose a directory
## you can set the directory name beforehand like this:
## >> basedir="../competitors/data"; dirname="C1"; evaalscore

## The necessary files are 'groundTruth.txt' and the file 'distance', which
## contains one line of text: either "localdistance" or "wgs84distance",
## depending on the reference system used for the coordinates.  It must
## be either a local system with unit in meters or else WGS84. In the
## latter case, the coordinates must be in the order long, lat.

## The files 'buttonsPressed.log' and 'positions.log' are created by the
## stepLogger app.  If they do not exist, and dopath=true is set when calling
## evaalscore, two dummy random files are created in their place which emulate
## an erro-free run.  The create_dummy_positions functions can be used to
## create a synthetic positions file with random errors for testing.

## If files relative to a defined scenario are present, the relative
## report is created.  Maps should be PNG 8 bit image files used as the
## graph backgrounds, '*.corners' should be text files containing the
## ref system coordinates of the four corners of the map in the usual
## image order: NO, NE, SE, SW.

################################################################
## TODO:
## arrow pointing outwards for big errors
## add more checks
## set titles of pdf files
## do not use plot for maps and arrows: draw everythng by hand
## movie: starting title slide
## movie: trim borders
## movie: dynamic stats computation
## movie: clock sound at every frame
## movie: speed up by removing repeated plot
## movie: writeVideo from pkg video does not work with svg

## Scenarios are tried in order starting from the first one.

scenarios = struct();

################ Empty scenario
s = struct();
s.reporttitle = "Empty scenario";
s.floorpenalty = 15;	# Floor penalty [m] used when computing distance

## Maps info
s.paperorientation = 'landscape';
s.paperposition = [0.5 0.5 29 21];	# centimetres
s.spsr = 2;				# subplot rows
s.spsc = 3;				# subplot columns
s.spi_legax = 0;			# subplot index of legend axes (0 means use file)
s.spi_legend = s.spsc;			# subplot index for legend and stats
s.spp_legend = [];			# explicit legend position
s.spi_distrib = s.spsr*s.spsc;		# subplot index for distribution graph
s.spp_distrib = [];			# distribution graph position
s.statx = 0;				# statistics xpos relative to distr. graph
s.staty = 0;				# statistics ypos relative to distr. graph
s.mapinfo = struct("map", {"F1", "F0"	# file basename
			  },
		   "rm", {0, 0},	# no scaling
		   "title", {"Empty first floor", "Empty ground floor"},
		   "floor", {1, 0},	# floor number
		   "spi", {1:2, 4:5},	# subplot index
		   "spp", {[],[]},	# subplot explicit position
		   ## Do not touch these
		   "psz", [],		# pixel width and height
		   "p2r", [],		# pixel to ref tform
		   "ah", []		# axes handle
		  );

scenarios.empty = s; clear s;

################ IPIN 2023 scenario (show only floors 0 and -2)
s = struct();
s.reporttitle = "IPIN competition 2023, Nürnberg (DE)";
s.floorpenalty = 15;	# Floor penalty [m] used when computing distance

## Maps info
s.paperorientation = 'landscape';
s.paperposition = [1 1.5 27 19];	# centimetres
s.spsr = 2;				# subplot rows
s.spsc = 3;				# subplot columns
s.spi_legax = 0;			# subplot index of legend axes (0 means use file)
s.spi_legend = s.spsc;			# subplot index for legend and stats
s.spp_legend = [0.75 0.75 0.20 0.20];	# explicit legend position
s.spi_distrib = s.spsr*s.spsc;		# subplot index for distribution graph
s.spp_distrib = [0.75 0.05 0.20 0.20];	# distribution graph position
s.statx = -4;				# statistics xpos relative to distr. graph
s.staty = 2.8;				# statistics ypos relative to distr. graph
s.mapinfo = struct("map", {"P0_MIC", "P2_MIC" # file basename
			  },
		   "rm", {0, 0},	# index of reference map for scaling
		   "title", {"Ground floor", "Underground floor (-2)"},
		   "floor", {0, -2},	# floor number
		   "spi", {1:2, 4:5},	# subplot index
		   "spp", {[.05 .45 .60 .60], [.05 .00 .60 .60]},
		   "psz", [], "p2r", [], "ah", [] # internal use
		  );

scenarios.evaal2023_2floors = s; clear s;

################ IPIN 2019 scenario
s = struct();
s.reporttitle = "IPIN competition 2019, Pisa (IT)";
s.floorpenalty = 15;	# Floor penalty [m] used when computing distance

## Maps info
s.paperorientation = 'landscape';
s.paperposition = [0.3 0 30 21];	# centimetres
s.spsr = 3;				# subplot rows
s.spsc = 3;				# subplot columns
s.spi_legax = 0;			# subplot index of legend axes (0 means use file)
s.spi_legend = 1;			# subplot index for legend and stats
s.spp_legend = [];			# explicit legend position
s.spi_distrib = 3;			# subplot index for distribution graph
s.spp_distrib = [];			# distribution graph position
s.statx = -30;				# statistics xpos relative to distr. graph
s.staty = 1;				# statistics ypos relative to distr. graph
s.mapinfo = struct("map", {"F2-CNR", "F1-CNR", "F0-CNR" # file basename
			  },
		   "rm", {0, 0, 0},	# no scaling
		   "title", {"Second floor", "First floor", "Ground floor"},
		   "floor", {2, 1, 0},	# floor number
		   "spi", {[4 7], [5 8], [6 9]}, # subplot
		   "spp", {[.0 .1 .31 .5], [.33 .1 .31 .5], [.66 .1 .31 .5]},
		   "psz", [], "p2r", [], "ah", [] # internal use
		  );

scenarios.evaal2019 = s; clear s;

################ IPIN 2017 scenario
s = struct();
s.reporttitle = "IPIN competition 2017, Sapporo (JP)";
s.floorpenalty = 15;	# Floor penalty [m] used when computing distance

## Maps info
s.paperorientation = 'landscape';
s.paperposition = [1 1.5 27 19];	# centimetres
s.spsr = 2;				# subplot rows
s.spsc = 3;				# subplot columns
s.spi_legax = 0;			# subplot index of legend axes (0 means use file)
s.spi_legend = s.spsc;			# subplot index for legend and stats
s.spp_legend = [0.75 0.75 0.20 0.20];	# explicit legend position
s.spi_distrib = s.spsr*s.spsc;		# subplot index for distribution graph
s.spp_distrib = [0.75 0.05 0.20 0.20];	# distribution graph position
s.statx = -4;				# statistics xpos relative to distr. graph
s.staty = 2.8;				# statistics ypos relative to distr. graph
s.mapinfo = struct("map", {"P1-HOK", "P0-HOK" # file basename
			  },
		   "rm", {0, 1},	# index of reference map for scaling
		   "title", {"First floor", "Ground floor"},
		   "floor", {1, 0},	# floor number
		   "spi", {1:2, 4:5},	# subplot index
		   "spp", {[.15 .54 .45 .45], [.15 .00 .45 .45]},
		   "psz", [], "p2r", [], "ah", [] # internal use
		  );

scenarios.evaal2017 = s; clear s;

################ EvAAL-IPIN 2016 scenario
s = struct();
s.reporttitle = "EvAAL-IPIN competition 2016, Alcalá de Henares (ES)";
s.floorpenalty = 15;	# Floor penalty [m] used when computing distance

## Maps info
s.paperorientation = 'landscape';
s.paperposition = [1 1.5 27 19];	# centimetres
s.spsr = 2;				# subplot rows
s.spsc = 5;				# subplot columns
s.spi_legax = 0;			# subplot index of legend axes (0 means use file)
s.spi_legend = s.spsc;			# subplot index for legend
s.spp_legend = [0.84 0.75 0.20 0.20];	# explicit legend position
s.spi_distrib = s.spsr*s.spsc;		# subplot index for distribution graph
s.spp_distrib = [0.84 0.05 0.20 0.20];	# distribution graph position
s.statx = -4;				# statistics xpos relative to distr. graph
s.staty = 2.8;				# statistics ypos relative to distr. graph
s.mapinfo = struct("map", { "P3_UAH",  "P2_UAH" # file basename
			    "P0_UAH",  "P1_UAH"
			  },
		   "rm", {0, 0; 0, 0},	# no map scaling
		   "title", {"Third floor",	"Second floor"
			     "Ground floor",	"First floor"
			    },
		   "floor", {3, 2; 0, 1}, # floor number
		   "spi", {1:2, 3:4	# subplot index
			   6:7, 8:9
			  },
		   "spp", {		# subplot explicit position
			   [.00 .50 .37 .37], [.42 .50 .37 .37],
			   [.00 .00 .37 .37], [.42 .00 .37 .37]},
		   "psz", [], "p2r", [], "ah", [] # internal use
		  );

scenarios.evaal2016 = s; clear s;

################ EvAAL-ETRI 2015 scenario
s = struct();
s.reporttitle = "EvAAL-ETRI competition 2015, Banff (CA)";
s.floorpenalty = 15;	# Floor penalty [m] used when computing distance

## Maps info
s.paperorientation = 'portrait';
s.paperposition = [];			 # centimetres
s.spsr = 3;				# subplot rows
s.spsc = 4;				# subplot columns
s.spi_legax = 0;			# subplot index of legend axes (0 means use file)
s.spi_legend = s.spsc;			# subplot index for legend
s.spp_legend = [0.75 0.75 0.20 0.20];	# explicit legend position
s.spi_distrib = s.spsr*s.spsc;		# subplot index for distribution graph
s.spp_distrib = [0.77 0.05 0.20 0.20];	# distribution graph position
s.statx = 2;				# statistics xpos relative to distr. graph
s.staty = 2;				# statistics ypos relative to distr. graph
s.mapinfo = struct("map", { "KCCI-L3",  "PDC-L3" # file basename
			    "KCCI-L2",  "PDC-L2"
			    "KCCI-L1",  "PDC-L1"
			  },
		   "rm", {0, 0; 0, 0; 0, 0},	# no map scaling
		   "title", {"maps not to scale",	""
			     "",			""
			     "",			""
			    },
		   "floor", {3, 3; 2, 2; 1, 1},	# floor number
		   "spi", {1:2, 3		# subplot index
			   5:6, 7
			   9:10, 11
			  },
		   "spp", {		# subplot explicit position
			   [.00 .64 .45 .29], [.45 .64 .30 .29],
			   [.00 .32 .45 .29], [.45 .32 .30 .29],
			   [.00 .00 .45 .29], [.45 .00 .30 .29]},
		   "psz", [], "p2r", [], "ah", [] # internal use
		  );

scenarios.evaal2015 = s; clear s;

################ End of scenarios

#

## Work around some Octave imread bugs
function img = fileimread(file)
  img = imread(file);
  if (size(img,3) == 3 && islogical(img))	# work around Octave bug #49140
    img = double(img);				# convert to double
#  elseif (isgray(img))				# work around Octave bug #49137
#    img = repmat(img, 1, 1, 3);		# convert to RGB
  endif
endfunction


## Draw the path only, without errors
function draw_path (scnr, directory)
  gtruthfile = [directory "/groundTruth.txt"];
  if !exist(gtruthfile, 'file')
    error("%s does not exist: cannot draw path\n", gtruthfile);
  endif

  buttonsfile = [directory "/buttonsPressed.log"];
  positionsfile = [directory "/positions.log"];
  pathfile = [directory "/path.pdf"];
  buttonsfilesaved = [buttonsfile ".saved"];
  positionsfilesaved = [positionsfile ".saved"];

  for savedfn = {buttonsfilesaved, positionsfilesaved}
    sfn = savedfn{:};
    if exist(sfn, 'file')
      error("found file '%s': check it, remove it and try again\n", sfn);
    endif
  endfor

  unwind_protect
    bfren = pfren = false;
    bfexist = exist(buttonsfile, 'file');
    pfexist = exist(positionsfile, 'file');
    if (bfexist) rename(buttonsfile, buttonsfilesaved); bfren = true; endif
    if (pfexist) rename(positionsfile, positionsfilesaved); pfren = true; endif

    ## Create a dummy buttonsPressed.log
    create_dummy_buttonpresses(directory, buttonsfile);

    ## Create path-only output
    create_dummy_positions(directory, 0);
    process_scenario(scnr, directory, true, false);
    print(pathfile);

  unwind_protect_cleanup
    if (bfren) movefile(buttonsfilesaved, buttonsfile, 'f') # overwrite destination
    else warning("creating dummy '%s' as one is missing\n", buttonsfile)
    endif
    if (pfren) movefile(positionsfilesaved, positionsfile, 'f') # overwrite destination
    else warning("creating dummy '%s' as one is missing\n", positionsfile)
    endif
  end_unwind_protect

endfunction


## Create a random "buttonsPressed.log" file useful for testing
##
## Call like this:
##  create_dummy_buttonpresses ("testevaal")
function  create_dummy_buttonpresses (directory, buttonsfile)
  if (nargin == 1)
    buttonsfile = [directory "/buttonsPressed.log"];
  endif
  gtruthfile = [directory "/groundTruth.txt"];

  gth = fopen(gtruthfile, 'r');
  out = textscan(gth, "%s%f%f%f");
  fclose(gth);
  [l x y f] = out{:};
  bfh = fopen(buttonsfile, 'w');
  for t = 1:length(l)
    fprintf(bfh, "%ld : %s : 0 : 0 : 0\n", 1000*(fix(time())+10*t), l{t});
  endfor
  fclose(bfh);
endfunction


## Create a random "positions.log" file useful for testing
## N is how many new pos to create between existing ones.
## If N is 0, create perfect positions, that is, without errors
##
## Call like this:
##  create_dummy_positions ("testevaal", 5)
function  create_dummy_positions (directory, n)
  positionsfile = [directory "/positions.log"];

  t = textread([directory "/buttonsPressed.log"], "%f%*%*%*%*", 'delimiter', ":");
  [l x y f] = textread([directory "/groundTruth.txt"], "%s%f%f%f");
  if (length(t) != length(l))
    warning("%d ground truth labels, but %d button taps\n",
	  length(l), length(t));
    return
  endif

  if (n == 0)			# copy reference points
    spd = [t-1 x y f];
  else				# interpolate and perturbate ref points
    p = interp1([t x y f], 0:1/(n+1):length(t)+2, 'extrap');
    absdiff = quantile(abs(diff(p)), 0.5)
    d = n * (rand(size(p))-1) .* absdiff;
    pd = p+d;
    pd(:,4) = round(pd(:,4));		# floor must be integer
    spd = sortrows(pd);			# sort by increasing times
  endif

  fid = fopen(positionsfile, 'w');
  fprintf(fid, "%.0f\t%.6f\t%.6f\t%d\n", spd');
  fclose(fid);
endfunction

## Create the contents of an IMAGE.corners file from IMAGE.png and IMAGE.points
function mapcorners = create_corners(image, projection = 'affine', opt=0)
  pkg load image

  imagefile = [ image ".png"];
  pointsfile = [ image ".points"];
  cornersfile = [ image ".corners"];

  ## Create coordinate conversion matrix pixel_to_ref from pointsfile
  pw = imfinfo(imagefile).Width;
  ph = imfinfo(imagefile).Height;
  corners_pixel = [0,0; pw,0; pw,-ph; 0,-ph];
  [mapX mapY pixelX pixelY] = textread(pointsfile, "%f%f%f%f%*", 'delimiter', ",", 'headerlines', 1);
  pixel_to_ref = cp2tform([pixelX pixelY], [mapX mapY], projection, opt);

  ## Create image corners cordinates and store them into cornersfile
  mapcorners = tformfwd(pixel_to_ref, corners_pixel);
  if (nargout == 0)
    printf("Contents to be saved into %s:\n", cornersfile);
    printf("%12.6f%12.6f\n", mapcorners');
  endif

  if exist(cornersfile, 'file')
    error("File %s exists: contents not overwritten\n", cornersfile);
  else
    save_precision(10, 'local');
    save("-text", cornersfile, "mapcorners");
  endif

endfunction

#

## Look in BASEDIR for directories containing the necessary files and
## ask the user to choose one of them, return the chosen dir as a string,
## that is the dir name relative to BASEDIR
function dir = chosendir(basedir);

  ## The files expected in a dir with one or no maps.  Optional files "map"
  ## and "map_corners" are required to add a background map
  files_single_map = { "buttonsPressed.log",
		       "positions.log",
		       "distance"};

  files = readdir(basedir);
  gooddiridx = [];
  for i = 1:length(files)
    if (isfolder([basedir "/" files{i}])
	&& !(strcmp(".", files{i}) || strcmp("..", files{i})))
      dirfiles = readdir([basedir "/" files{i}]);
      found_single = 0;
      for j = 1:length(dirfiles)
	found_single += any(strcmp(dirfiles{j}, files_single_map));
      endfor
      if (found_single == length(files_single_map))
	gooddiridx(end+1) = i;
      endif
    endif
  endfor
  dirs = files(gooddiridx);

  assert(length(dirs) > 0, "no good directories found under '%s', exiting\n", basedir);
  choice = 0;
  choice = menu("Choose a directory", {"no choice: EXIT" dirs{:}});
  if (choice < 2)
    dir = ""
  else
    dir = dirs{choice-1}	# first choice is "no choice"
  endif

endfunction


## Convert vectors of X and Y coordinates to pixel coords
function [px py] = to_pixel_coords (ph, pixel_to_ref, x, y)
  p = tforminv(pixel_to_ref, [x y]);
  px = p(:,1);
  # As usual, we set the image origin in the top-left corner but, like
  # QGIS' georeferencer, we do not follow the usual convention of y axis
  # going downwards, consequently y coordinates are negative.
  # To adjust things, below we use flipud on the image.
  py = p(:,2) + ph;
endfunction


## Distance in meters given longitude and latitude
## Points P1 and P2 are [long lat floor] or nx3 matrices
function d = wgs84distance(p1, p2, floorpenalty = 0)
  lng1 = p1(:,1); lat1 = p1(:,2);
  lng2 = p2(:,1); lat2 = p2(:,2);
  dLng = lng2 - lng1;
  dLat = lat2 - lat1;

  ## Equirectangular approximation
  c = norm(pi/180*[dLng .* cosd(lat1), dLat], 'rows');
  d = 6371000 * c;			# 6371000 is the mean Earth radius [m]

  ## Add floor penalty
  d += floorpenalty * abs(p1(:,3) - p2(:,3));
endfunction


## Distance given coordinates in their measurement unit
## Points P1 and P2 are [x y f] or nx3 matrices
function d = localdistance(p1, p2, floorpenalty = 0)
  xydiff = (p2-p1)(:,1:2);		# diffs of xy coordinates
  d = norm(xydiff, 'rows');

  ## Add floor penalty
  d += floorpenalty * abs(p1(:,3) - p2(:,3));
endfunction

#

## Compute path statistics
function create_path_stats (directory, distance)
  gtruthfile = [directory "/groundTruth.txt"];
  if !exist(gtruthfile, 'file')
    warning("%s does not exist: cannot create path statistics\n", gtruthfile);
  else
    [l x y f] = textread(gtruthfile, "%s%f%f%f");
    floorstart = 1 + [0; find(diff(f))];
    floorend = [find(diff(f)); length(f)];
    sumfloors = 0;
    sumpoints = 0;
    sumlen = 0;
    printf("\nFloor   | horizontal path stats | floor transition \n%s", repmat("-",1,52));
    for fidx = 1:length(floorstart);
      sidx = floorstart(fidx);
      eidx = floorend(fidx);
      spos = [x y f](sidx:eidx-1, :);
      epos = [x y f](sidx+1:eidx, :);
      points = eidx - sidx + 1;
      #wgs84distance(spos, epos)
      len = sum(feval(distance, spos, epos));
      printf("\n%3d     |%4d points |%7.1f m",
	     f(floorstart(fidx)), points, len );
      sumfloors += 1;
      sumpoints += points;
      sumlen += len;
      if (eidx < length(x))	# if floor change ahead
	len = feval(distance, [x y f](eidx, :), [x y f](eidx+1, :));
	printf(" |%3d -> %2d  |%5.1f m",
	       f(eidx), f(eidx+1), len );
	sumlen += len;
      endif
    endfor
    printf("\nTotals  -------------------------\n%d paths |%4d points |%7.1f m\n",
	   sumfloors, sumpoints, sumlen);
  endif

  buttonsfile = [directory "/buttonsPressed.log"];
  if exist(buttonsfile, 'file')
    [bpt bpl bpx bpy bpf] = textread(buttonsfile, "%f%s%f%f%f", 'delimiter', ":");
    sumtime = ceil(range(bpt) / 1000);
    printf("\nTotal time taken %d'%d\"",
	   fix(sumtime/60), rem(sumtime, 60));
    if exist("sumlen", 'var')
      printf(", average speed %.1f m/s", sumlen/sumtime);
    endif
    printf("\n");
  endif

  printf("\n");

endfunction


#

## Compute point errors, possibly draw scenario graph or make movie
## Read the data files, put the reference coordinates in RX RY and the
## estimated coordinates in EX EY, draw the estimated path with the map image
## in background, if available
function [rl rt rx ry rf et ex ey ef] = process_scenario (s, directory, doplot, domovie, domovie_using_mkmovie)

  buttonsfile = [directory "/buttonsPressed.log"];
  positionsfile = [directory "/positions.log"];
  groundtruthfile = [directory "/groundTruth.txt"];

  ## Read button pressed time, button label, reference coordinates (ground truth)
  [bpt bpl bpx bpy bpf] = textread(buttonsfile, "%f%s%f%f%f", 'delimiter', ":");

  assert(issorted(bpt), "non-increasing times in %s", buttonsfile);

  ## If the ground truth file is present, the bpx and bpy coordinates that
  ## we have just read are dummy ones, and we should recover them from the
  ## ground truth file; also check that the dummy values are 0
  assert(all([bpx; bpy; bpf] == 0),
	 "%s is present, but not all coordinates in %s are 0\n",
	 groundtruthfile, buttonsfile);
  [gtl gtx gty gtf] = textread(groundtruthfile, "%s%f%f%f");

  ## Check labels and floor
  gtl = strtrim(gtl); bpl = strtrim(bpl); # clean up labels
  assert(numel(gtl) == numel(unique(gtl)), "non-unique labels in %s", groundtruthfile);
  assert(numel(bpl) == numel(unique(bpl)), "non-unique labels in %s", buttonsfile);
  assert(all(round(gtf) == gtf), "non-integer floor in %s", groundtruthfile);

  ## Find indices of button pressures into ground truth by matching
  ## the button pressure labels bpl to ground truth labels gtl
  [found idx] = ismember(bpl, gtl);
  assert(all(found), "%s does not contain all labels in %s\n",
	 groundtruthfile, buttonsfile);
  assert(issorted(idx), "labels in %s not sorted as %s", buttonsfile, groundtruthfile);
  if (length(bpl) != length(gtl))
    warning("%d ground truth labels, but %d button pressures\n",
	    length(gtl), length(bpl));
  endif

  ## Set reference coordinates of button pressures
  bpx = gtx(idx); bpy = gty(idx); bpf = gtf(idx); # button pressures coordinates
  bpl = gtl(idx);				  # button pressure labels

  rl = bpl; rt = bpt;			# OUTPUT labels and time
  rx = bpx; ry = bpy; rf = bpf;		# OUTPUT reference coordinates

  ## Read estimation time, estimated coordinates (competitor's output)
  [pt px py pf] = textread(positionsfile, "%f%f%f%f");

  ## Check floor
  assert (all(round(pf) == pf), "non-integer floor in %s", positionsfile);

  ## Check for swapped lat<->long (should be longitude,latitude)
  mgt = mean([gtx gty]);
  assert(norm([px py]-mgt) < norm([py px]-mgt), "swapped longitude and latitude\n");

  ## Check that times are (almost) increasing and possibly correct them
  maxdecrease = 550;			# milliseconds
  minptdiff = min(diff(pt));
  if (minptdiff < 0)
    if (-minptdiff < maxdecrease)
      warning("non-increasing times in %s: max decrease %d ms\n",
	      positionsfile, -minptdiff);
      [pt idx] = sort(pt);
      px = px(idx); by = py(idx); pf = pf(idx);
    else
      error("non-increasing times in %s: max decrease %d > %d ms\n",
	    positionsfile, -minptdiff, maxdecrease);
    endif
  endif

  ## Check that first pos precedes first button pressure,
  ## possibly correct adding one pos at the beginning
  if (pt(1) > bpt(1))
    warning("%s starts after first button press, adding one initial log position\n",
	    positionsfile);
    pt = [bpt(1); pt];
    px = [px(1); px]; py = [py(1); py]; pf = [pf(1); pf];
  endif

  pidx = lookup(pt, bpt);			  # indices of competitor estimates
  bex = px(pidx); bey = py(pidx); bef = pf(pidx); # estimates at time of button pressed
  et = pt(pidx); ex = bex; ey = bey; ef = bef;	  # OUTPUT estimated coordinates
  assert(all(rt >= et), "internal consistency check failed");

  if (!doplot) return; endif;	# stop here after reading data files


  ## Set up the map images
  #########################################################################
  alpha = 0.7;				# transparency of map image
  light = 0.9;				# light level of grey behind transparent map image

  mapno = numel(s.mapinfo);		# number of maps
  for m = 1:mapno			# for all the maps

    ## Create coordinate conversion matrix pixel_to_ref
    mapfile = [directory "/" s.mapinfo(m).map ".png"];
    cornersfile = [directory "/" s.mapinfo(m).map ".corners"];
    pw = imfinfo(mapfile).Width;
    ph = imfinfo(mapfile).Height;
    s.mapinfo(m).psz = [pw ph];
    corners_pixel = [0,0; pw,0; pw,-ph; 0,-ph];
    corners_ref = load('-ascii', cornersfile); # Ignore Octave headers
    pixel_to_ref = cp2tform(corners_pixel, corners_ref, 'projective');
    s.mapinfo(m).p2r = pixel_to_ref; # store it in s.mapinfo

    ## Draw the map as background
    mapimg = alpha*light*255 + (1-alpha)*fileimread(mapfile);
    ah = subplot(s.spsr, s.spsc, s.mapinfo(m).spi);
    if !isempty(s.mapinfo(m).spp)	# spp contains explicit position info
      ## Scale and center the map with respect to another map used as reference
      rm = s.mapinfo(m).rm;		# index of reference map
      if (s.mapinfo(m).rm != 0)	# scale the map
	scale = s.mapinfo(m).psz ./ s.mapinfo(rm).psz; # [w h] scale
	s.mapinfo(m).spp(3:4) .*= scale;
	s.mapinfo(m).spp(1) -= s.mapinfo(m).spp(3) * (scale(1)-1)/2;
      endif
      set(ah, 'position', s.mapinfo(m).spp);
    endif
    s.mapinfo(m).ah = ah;			# axis handle

    ## We need axis('xy') to be able to plot later on the same graph,
    ## that's why we use flipud and 'xy' rather than simply axis('ij')
    imshow(flipud(mapimg));
    axis('xy', 'image', 'off', 'nolabel');
    hold on
  endfor

  ## Draw the estimated path, estimated points, real points, error segments
  #########################################################################

  ## Create coordinates of error lines for plotting, NA is used to create
  ## separate segments when plotting
  errx = reshape([bpx ex NA(size(bpx))]', [], 1);
  erry = reshape([bpy ey NA(size(bpy))]', [], 1);
  wfidx = (bpf != bef);			# indices of bad floor estimates

  ## Do the plotting, for each img do the plot for all points in the
  ## corresponding floor
  printf("Creating %s:", {"plot", "movie"}{domovie+1});

  ## The number of points is N = length(px), the number of maps is mapno
  ## (N,mapno) is the size of pp[xy], pbr[xy], pbe[xy], pep[xy]
  ## (3*N,mapno) is the size of perr[xy]
  ## (N-1,mapno) is the size of dep[xy], dlen
  ## (N-1,1) is the size of depf
  for m = 1:mapno
    ph = s.mapinfo(m).psz(2);		# map height in pixels
    pixel_to_ref = s.mapinfo(m).p2r;

    [ppx(:,m), ppy(:,m)] = to_pixel_coords(ph, pixel_to_ref, px, py);
    [perrx(:,m), perry(:,m)] = to_pixel_coords(ph, pixel_to_ref, errx, erry);
    [pbrx(:,m), pbry(:,m)] = to_pixel_coords(ph, pixel_to_ref, bpx, bpy);
    [pbex(:,m), pbey(:,m)] = to_pixel_coords(ph, pixel_to_ref, bex, bey);
  endfor

  pepx = pbex;		 # pbex is estimated path, use ppx for whole path
  pepy = pbey;		 # pbey is estimated path, use ppy for whole path
  depx = diff(pepx); depy = diff(pepy); depf = diff(bpf);
  dlen = sqrt(sumsq(cat(3,depx,depy),3)); # used to set the arrow head sizes below

  msize = 2.5;				# size of markers
  lwidth = 1;				# width of thick lines
  hsize = 6*lwidth;			# size of arrow heads
  if strcmp(graphics_toolkit(), 'gnuplot')
    lwidth = 4;
  endif

  hold on

  for m = 1:mapno
    fi = (rf == s.mapinfo(m).floor);	# index of floor's points

    ## Ground truth
    plot(s.mapinfo(m).ah,
    	 pbrx(fi,m), pbry(fi,m), "og", 'markersize', msize, 'markerfacecolor', 'green');
    axis('equal');
  endfor

  if (domovie)
    set(1,				# increase figure resolution for PNG
	'papertype', '<custom>', 'paperorientation', 'landscape',
	'paperunits', 'points', 'paperposition', [0 0 1920 1080]);

    if domovie_using_mkmovie
      moviefn = fullfile(directory, "path.ogg");
      mkmovie('init', moviefn);
    else
      moviefh = VideoWriter(fullfile(directory, "path.mp4"));
      moviefh.FrameRate = 1.2;
      open (moviefh);
    endif
    for f = 1:rows(dlen)		# frames
      pf = f-1;				# previous frame
      chr = 0;				# real circle handle
      for m = 1:mapno			# maps
	if (s.mapinfo(m).floor == rf(f)) # point f is on same floor as map m
	  mcol = 'blue';
	  if wfidx(f)
	    mcol = 'red';
	  endif
	  ## TODO: no need to call plot repeatedly, just do it once and change elements
	  ## See https://stackoverflow.com/questions/54870121/how-to-create-a-mp4-video-out-of-octave-via-plot
	  plot(s.mapinfo(m).ah,
	       perrx(3*f-2+(0:1),m), perry(3*f-2+(0:1),m), "-r",
	       pbex(f,m), pbey(f,m), "o",
	       'markeredgecolor', mcol, 'markerfacecolor', mcol, 'markersize', msize);

	  ## Using rectangle() with total curvature to draw circles
	  r = min(s.mapinfo(m).psz)/15;	  # radius
	  xc = pbrx(f,m); yc = pbry(f,m); # centre
	  chr = rectangle('position', [xc-r yc-r 2*r 2*r], 'curvature', [1 1],
			  'edgecolor', 'none', 'facecolor', 'green', 'facealpha', 0.3);
	  r /= 2;			  # radius
	  xc = pbex(f,m); yc = pbey(f,m); # centre
	  che = rectangle('position', [xc-r yc-r 2*r 2*r], 'curvature', [1 1],
			  'linewidth', 3, 'edgecolor', 'red', 'edgealpha', 0.5);

	endif
	if (pf > 0 && dlen(pf,m) > 0
	    && s.mapinfo(m).floor == rf(pf)) # point pf is on same floor as map m
	  arrowcolour = 'blue';
	  if depf(pf)
	    arrowcolour = 'cyan';
	  endif
	  quiver(pepx(pf,m), pepy(pf,m), depx(pf,m), depy(pf,m), 0, 'color', arrowcolour,
  		 'linewidth', lwidth, 'maxheadsize', hsize./dlen(pf));
	endif
      endfor
      legend('off');
      c = ".";  if (rem(f,10) == 0) c = "|"; endif
      printf(c);
      if domovie_using_mkmovie
	mkmovie('add', moviefn, pf);	# add this frame to movie
      else
	## drawnow; writeVideo(moviefh, getframe(gcf)); # faster, but does not work with gnuplot
	## Cannot use getframe with gnuplot
	tmpimgfn = "tmpimage.png";
	saveas(1, tmpimgfn);
	writeVideo(moviefh, imread(tmpimgfn));
      endif
      if (chr != 0)			# if circles have been drawn
	delete(chr); delete(che);	# delete them
      endif
    endfor
    ## TODO: put these inside an unwind_protect_cleanup
    if domovie_using_mkmovie
      mkmovie('close', moviefn, 1.2);
    else
      close(moviefh);
      delete(tmpimgfn);
    endif

  else
    set(1,				# landscape on A4 paper
	'papertype', 'a4', 'paperorientation', s.paperorientation,
	'paperunits', 'centimeters', 'paperposition', s.paperposition);

    for m = 1:mapno
      printf([" " s.mapinfo(m).map]);	# show advancement
      fflush(1);			# necessary for Windows after printf

      fi = (rf == s.mapinfo(m).floor);	# index of floor's points
      fi3 = repmat(fi', 3, 1);		# repeat each element three times
      px = pepx(:,m); py = pepy(:,m);
      px(!fi) = py(!fi) = NA;

      plot(s.mapinfo(m).ah,
	   px, py, "-b;estimated path;", 'linewidth', lwidth,
	   perrx(fi3,m), perry(fi3,m), "-r;errors;",
	   pbex(fi,m), pbey(fi,m), "ob;estimated points;", 'markersize', msize,
	   pbex(fi&wfidx,m), pbey(fi&wfidx,m),
	   "or;wrong floor;", 'markerfacecolor', 'red', 'markersize', msize);

      ## Overwrite the estimated path plotted above with arrows.
      ## The head size is set for each arrow using the "hsize" variable.
      ## The plot instruction above must nonetheless include
      ## the estimated path so that it appears in the legend.
      for i = 1:rows(dlen)
	if (fi(i) && dlen(i,m) > 0)
	  arrowcolour = 'blue';		# path to same floor
	  if depf(i)
	    arrowcolour = 'cyan';	# path to a different floor
	  endif
	  quiver(pepx(i,m), pepy(i,m), depx(i,m), depy(i,m), 0, 'color', arrowcolour,
  		 'linewidth', lwidth, 'maxheadsize', hsize./dlen(i));
	endif
      endfor

      ## Plot accessory elements
      if isaxes(s.mapinfo(m).ah) && !isempty(s.mapinfo(m).title)
	fsiz = get(s.mapinfo(m).ah, 'fontsize') * 1.25;
	text(s.mapinfo(m).ah, mean(xlim()), mean(ylim())-[s.mapinfo.psz](2)/2-4*fsiz,
	     s.mapinfo(m).title,
	     'fontsize', fsiz, 'horizontalalignment', 'center');
      endif
      if (m == s.spi_legax)
	hleg = legend('location', 'northeastoutside');
	legpos = get(hleg, 'position');
	set(hleg, 'position', legpos+[0 -0.15 0 0]);
      else
	legend off;
      endif
    endfor

  endif				# movie or plot

  hold off
  printf(" done.\n");

endfunction


#
################################################################
## Main program
################################################################

## Set default values
exist("basedir") || (basedir = ".");
exist("doplot") || (doplot = true);
exist("dopath") || (dopath = false);
exist("domovie") || (domovie = false);
if domovie
  ## 2023-09 ffmpeg shipped with Octave on Windows cannot read svg
  domovie_using_mkmovie = (!ispc() && exist("mkmovie", 'file') == 2);
  if !domovie_using_mkmovie
    if isempty(pkg('list', 'video'))
      warning("'mkmovie.m' not present and package 'video' not installed.  Cannot make a movie.");
      domovie = false;
    else
      pkg load video
    endif
  endif
endif
(domovie) && (doplot = true);
if (! exist("dirname"))
  dirname = chosendir(basedir);
endif
assert(!isempty(dirname), "no directory chosen, exiting\n");

directory = fullfile(basedir, dirname);

distancefile = fullfile(directory, "distance");
legendfile = fullfile(directory, "legend.png");
errorfile = fullfile(directory, "error.report");
qgisfile = fullfile(directory, "qgis_errors.txt");
buttonsfile = fullfile(directory, "/buttonsPressed.log");
positionsfile = fullfile(directory, "/positions.log");

## Get the distance type
distance = textread(distancefile, "%s"){1};

## Print statistics
create_path_stats(directory, distance);

## Select scenario
for [s name] = scenarios
  notfound = false;
  for mapfile = {s.mapinfo.map}
    if !exist(fullfile(directory, [mapfile{:} ".png"]), 'file')
      notfound = true;
      break;
    endif
  endfor
  if (notfound)
    continue;					# check next scenario
  else
    printf(["Using scenario '" name "\"\n"]);
    scnr = s;
    break
  endif
endfor
if (notfound)
  printf(["Using no scenario\n"]);
  scnr = struct("floorpenalty", 0, "reporttitle", "");
  if (doplot)
    printf("No known scenario in the selected directory: setting doplot=false\n");
    doplot = false;
  endif
endif
fflush(1);			# necessary for Windows

if isempty(pkg('list', 'image'))
  warning("cannot plot: package 'image' not installed.  Generating report only.\n");
  doplot = false
endif

if (doplot)
  pkg load image				# for tform functions
  warning('off', 'Octave:GraphicsMagic-Quantum-Depth');

  ## 2023-09 libgl2ps cannot handle PNG when generating SVG
  if verLessThan('Octave', "9")
    warning('off', "Octave:gnuplot-graphics");
    graphics_toolkit gnuplot
  endif

  ## Prevent plotting, do that at the end only
  figure(1, 'visible', 'off');
  clf(1);

  ## Generate path file

  ## If either buttonsfile or positionsfile do not exist,
  ## call draw_path so to create dummy ones
  if !exist(buttonsfile, 'file') || !exist(positionsfile, 'file')
    dopath = true
  endif
  if (dopath) draw_path(scnr, directory); endif
endif

## Read data, possibly draw the path with background maps or make movie,
## return the reference coordinates and the estimated coordinates
[rl rt rx ry rf et ex ey ef] = process_scenario(scnr, directory, doplot, domovie, domovie_using_mkmovie);

## Compute error statistics
poserr = feval(distance, [rx ry rf], [ex ey ef], scnr.floorpenalty);
ninetiethperc = quantile(poserr, 90/100, 1, 7)	# ninetieth percentile of error
thirdq = quantile(poserr, 3/4, 1, 7)		# third quartile of error

## Generate score file
if (doplot && !domovie)

  lah = subplot(scnr.spsr, scnr.spsc, scnr.spi_legend);
  axis(lah, 'off');

  ## Legend from file
  if (exist(legendfile, 'file') && scnr.spi_legax == 0)
    isempty(scnr.spp_legend) || set(lah, 'position', scnr.spp_legend);
    imshow(fileimread(legendfile));
  endif

  ## Global title
  gth = annotation('textbox', [0.4 0.925 0.2 0.2],
		   'horizontalalignment', 'center', 'fontsize', 15,
		   'interpreter', 'none',
		   'string', dirname);
  text(0.5, 0.5, scnr.reporttitle)

  ## Distribution graph
  ah = subplot(scnr.spsr, scnr.spsc, scnr.spi_distrib);
  isempty(scnr.spp_distrib) || set(ah, 'position', scnr.spp_distrib);
  maxerr = 20;				# max error to plot
  errno = length(poserr);			# number of error
  stairs([0; sort(poserr)], (0:errno)/errno);
  title("distribution");
  axis([0 maxerr 0 1], 'square');
  yticks(0:0.25:1)

  grid on
  hold on
  plot(thirdq*[1 1], [0 1], 'linewidth', 2, 'color', 'red');

  ## Statistics
  xp = scnr.statx; yp = scnr.staty;		# x, y position
  il = 0.13;					# interline space
  ln = 0;					# line number
  text(xp, yp, "Error stats [m]", 'fontsize', 12);
  text(xp, yp-il*(ln+=2), sprintf("median: %.1f", median(poserr)));
  text(xp, yp-il*(ln+=1), sprintf("mean: %.1f", mean(poserr)));
  text(xp, yp-il*(ln+=1), sprintf("rms: %.1f", std(poserr, 1)));
  text(xp, yp-il*(ln+=2), sprintf("90^{th} perc.: %.1f", ninetiethperc));
  text(xp, yp-il*(ln+=2), sprintf("3^{rd} quartile: %.1f", thirdq),
       'fontsize', 12, 'fontweight', 'bold');

  ## Generate score file
  scorefile = fullfile(directory, "score.pdf");
  print(scorefile);

endif

## Print Qgis data
fid = fopen(qgisfile, 'w');
    fprintf(fid,"id;x;y;z;id1;x1;y1;z1\n");
for l = 1:length(rl)
    fprintf(fid,"%s;",rl{l});
    fprintf(fid,"%.6f;%6f;%d;", [rx ry rf](l,:));
    fprintf(fid,"%s",rl{l});
    if (rf(l) != ef(l)) fprintf(fid,"*") endif
    fprintf(fid,";%.6f;%6f;%d\n", [ex ey ef](l,:));
endfor
fclose(fid);

## Print error report
fid = fopen(errorfile, 'w');
fprintf(fid, [scnr.reporttitle "\n\n"]);
fprintf(fid, "\t\t\t================ %s ================\n\n", dirname);
fprintf(fid, "\t\t\t  === Third quartile of error: %.1f [m] ===\n\n", thirdq);
fprintf(fid, "Floor penalty is %.1f metres\n", scnr.floorpenalty);
fprintf(fid, "Error statistics [m]: median: %.1f  mean: %.1f  RMS: %.1f\n\n",
	median(poserr), mean(poserr), std(poserr, 1));
fprintf(fid, "\
Label\tReference coordinates\t\tEstimated coordinates\tError[m]  Error[s]\
\n\n");
for l = 1:length(rl)
  fprintf(fid, "%5s%14.6f%10.6f%3d\t%11.6f%10.6f%3d",
	  rl{l}, [rx ry rf ex ey ef](l,:));
  if (rf(l) != ef(l))
    fprintf(fid, "%+8.1f", poserr(l));
  else
    fprintf(fid, "%8.1f", poserr(l));
  endif
  fprintf(fid, "%10.1f\n", (rt(l)-et(l))/1000);
endfor
fclose(fid);

## Set showplot to false in production
showplot = false;
if (doplot)
  if (showplot)
    figure(1, 'visible', 'on');
  else
    close(1);
  endif
endif

## Cleanup
clear doplot dopath domovie basedir dirname


## Local Variables:
## fill-column: 110
## End:
