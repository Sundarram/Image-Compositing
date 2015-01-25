///
/// CS-583: Introduction to Computer Vision
/// https://www.cs.drexel.edu/~kon/introcompvis/
///
/// Project 1 - Homography
///
/// @author SUNDAR RAM SWAMINATHAN
/// @date DUE 10/22/2013 11:59PM
///


// Needed for reading files
import java.io.FileNotFoundException;
import java.util.*;
// Needed for the various Matrix and Vector operations
import no.uib.cipr.matrix.*;
// Needed for functions like Max and Min
import java.lang.Math.*;
// Needed for the .ini file parsing
import org.ini4j.Ini;




// ini file reading constants including filename and key names
static final String INI_FILE = "homography.ini";
static final String SOURCE_IMAGE = "sourceImage";
static final String SOURCE_POINTS = "sourcePoints";
static final String TARGET_POINTS = "targetPoints";
static final String TARGET_IMAGE = "targetImage";
static final String SOURCE_MASK = "sourceMask";
static final String OUTPUT_IMAGE = "outputFile";

// Homographic indexes
static final int X = 0;
static final int Y = 1;
static final int Z = 2;

// Whether or not to print debug information
static final boolean DEBUG = true;

// Size of red dot to draw. Used in showPoints
static final int DOT_SIZE = 10;

///
/// GLOBALS
///

// To keep track on what section-step we are
int sectionKey = 0;
int showingSection = 0;
int showingStep = 0;
int[] process; // 1 = composite, 2 = rectify
int totalSections;

// Globals for the name of the sections
String[] heading;
String[] sourceImagePath;
String[] sourcePointsPath;
String[] targetPointsPath;
String[] targetImagePath;
String[] savePath;

// Globals for the images to show
PImage[] sourceImage;
PImage[] targetImage;
PImage[] sourceMask;
PImage[] result;

Matrix[] sourcePoints;
Matrix[] targetPoints;


void setup() {
  size(400, 400); // size(screen.width, screen.height); // fullscreen
  noLoop(); // turn off re-drawing loop

  // Read .ini file
  Ini parameters = new Ini();
  try {
    BufferedReader reader = createReader(INI_FILE);
    if(reader == null) {
      throw new FileNotFoundException("ini file not found");
    }

    parameters.load(reader);

    totalSections = parameters.size();

    // allocate memory based on the number of sections
    heading = new String[totalSections];
    sourceImagePath = new String[totalSections];
    sourcePointsPath = new String[totalSections];
    targetPointsPath = new String[totalSections];
    targetImagePath = new String[totalSections];
    savePath = new String[totalSections];
    process = new int[totalSections];
    sourceImage = new PImage[totalSections];
    targetImage = new PImage[totalSections];
    sourceMask = new PImage[totalSections];
    result = new PImage[totalSections];
    sourcePoints = new Matrix[totalSections];
    targetPoints = new Matrix[totalSections];

    // Iterate through sections of file
    Set sections = parameters.entrySet();
    Iterator sectionIterator = sections.iterator();
    while(sectionIterator.hasNext()) {
      Map.Entry e = (Map.Entry)sectionIterator.next();
      Ini.Section section = (Ini.Section)e.getValue();

      // Read each section to get the various possible input values

      // Read section settings
      heading[sectionKey] = (String)e.getKey();
      sourceImagePath[sectionKey] = section.get(SOURCE_IMAGE);
      sourcePointsPath[sectionKey] = section.get(SOURCE_POINTS);
      targetPointsPath[sectionKey] = section.get(TARGET_POINTS);
      targetImagePath[sectionKey] = section.get(TARGET_IMAGE);
      savePath[sectionKey] = section.get(OUTPUT_IMAGE);

            
      sourceImage[sectionKey] = loadImage(sourceImagePath[sectionKey]);
      
      sourcePoints[sectionKey] = readPoints(sourcePointsPath[sectionKey]);
      targetPoints[sectionKey] = readPoints(targetPointsPath[sectionKey]);

      // Fork to do compositing or rectification based on which values are set in ini section

      if(targetImagePath[sectionKey] != null) {
        //Composite
        targetImage[sectionKey] = loadImage(targetImagePath[sectionKey]);
        process[sectionKey] = 1;
        result[sectionKey] = composite(sourceImagePath[sectionKey], sourcePointsPath[sectionKey],
                                     targetImagePath[sectionKey], targetPointsPath[sectionKey]);
                                       
                                       
                                      
      } else {
        // Rectify
        process[sectionKey] = 2;
        result[sectionKey] = rectify(sourceImagePath[sectionKey], sourcePointsPath[sectionKey],
                                     targetPointsPath[sectionKey]);
      }

      // Save image
      if(!savePath[sectionKey].equals("")) {
        result[sectionKey].save(savePath[sectionKey]);
        
      }

      println("Section completed");
      sectionKey++;
    }
  } catch(FileNotFoundException e) {
      System.err.println("Failure reading ini file: " + e.getMessage());
      exit();
  } catch(IOException e) {
    System.err.println("Failure reading ini file: " + e.getMessage());
    exit();
  }


  println("Done.");

}


void keyReleased() {
  // Move one step forward
  if ( (process[showingSection] == 1 && showingStep < 3) ||
       (process[showingSection] == 2 && showingStep < 1) ) {
         showingStep++;
  } else {
    showingStep = 0;
    showingSection++;
    if (showingSection == totalSections) {
      exit();
    }
  }
  
  redraw();
}

///
/// Draw function
///
void draw() {
   // Set the image to draw
  PImage i = sourceImage[showingSection];
  Matrix newPoints = new DenseMatrix(0,0);
  
  if ( process[showingSection] == 1 ) { // composite
    if ( showingStep == 0 ) {
      i = sourceImage[showingSection];
      newPoints = sourcePoints[showingSection].copy();
      println("[COMPOSITE] Now showing source image with points highlighted");
    } else if ( showingStep == 1 ) {
      i = targetImage[showingSection];
      newPoints = targetPoints[showingSection].copy();
      println("[COMPOSITE] Now showing target image with points highlighted");
    } else if ( showingStep == 2 ) {
      //i = sourceMask[showingSection];
      i = result[showingSection]; // change this
      println("[COMPOSITE] Now showing mask image with points highlighted");
    } else if ( showingStep == 2 | showingStep == 3 ) {
      i = result[showingSection];
      println("[COMPOSITE] Now showing result image");
    }
  } else { //retify
    if ( showingStep == 0 ) {
      i = sourceImage[showingSection];
      println("[RETIFY] Now showing source image");
    } else if ( showingStep == 1 ) {
      i = result[showingSection];
      println("[RETIFY] Now showing result image");
    }
  }
  
  background(0);

  float divideBy = 1.0;
  if(i.width > width || i.height > height) {
      divideBy = max(float(i.width)/width, float(i.height)/height);
      newPoints.scale(1.0/divideBy);
  }

  int originX = floor((width-i.width/divideBy)/2);
  int originY = floor((height-i.height/divideBy)/2);

  for(int r = 0; r < newPoints.numRows(); r++) {
    newPoints.set(r, X, newPoints.get(r, X) + originX);
        newPoints.set(r, Y, newPoints.get(r, Y) + originY);
  }

  // Show the picture scaled appropriately
  image(i, originX, originY, i.width/divideBy, i.height/divideBy);

  // Show the points
  showPoints(newPoints);
}

///
/// Given an input image and pairs of corresponding points in files, deforms the image
///
PImage rectify(String sourceImagePath, String sourcePointsPath, String targetPointsPath) {
  // Load input image, source points and target points
  PImage sourceImage = loadImage(sourceImagePath);
  Matrix sourcePoints = readPoints(sourcePointsPath);
  Matrix targetPoints = readPoints(targetPointsPath);

 // Compute homography matrix and its inverse

  Matrix H    = computeH(sourcePoints, targetPoints);
  println("Homography matrix");
  println(H);
  
  Matrix Hinv = computeH(targetPoints, sourcePoints);
  println("Inverse homography matrix");
  println(Hinv);
  
 // Warp the image and return it

PImage result = warpImage(sourceImage, H, Hinv);
  return result;
}


PImage warpImage(PImage sourceImage, Matrix H, Matrix H_i) {
  // Determine the size of the new image
  no.uib.cipr.matrix.Vector origin = new DenseVector(2);
  no.uib.cipr.matrix.Vector dims = new DenseVector(2);
  newImageSize(sourceImage, H, origin, dims);

  println("New origin: " + origin);
  println("New image size: " + dims);

  // Create new image
  PImage newImage = createImage((int)Math.ceil(dims.get(X)), (int)Math.ceil(dims.get(Y)), RGB);
  int w = (int) Math.ceil(dims.get(X));
  int h = (int) Math.ceil(dims.get(Y));
    
  // Loop through the target image and determine which pixels are defined in the source image
  // - In the loop apply the inverse homography to get the image coordinates in the source image
  //   (For this, implement and use 'applyHomography')
  // - If applyHomography returns image coordinates that are defined in the source image, get the
  //   pixel value in the source image of that pixel using bilinear interpolation
  //   (For this, implement and use getColorBilinear)

  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      no.uib.cipr.matrix.Vector p = new DenseVector(new double[] {x, y});
      p.add(origin);
      no.uib.cipr.matrix.Vector q = applyHomography(H_i, p);
      color c = getColorBilinear(sourceImage, q);
      newImage.set(x, y, c);
    }
  }
  return newImage;
}


///
///  Given a homography H and an image, computes the bounding box of the image after applying H.
///
///  @param image - Source image
///  @param H - 3x3 homography matrix
///  @param origin - 2D vector defining upper left of new image
///  @param dimensions - 2D vector defining size of new image
///  @post index 0 and 1 of origin and dimensions have been set
///
void newImageSize(PImage image, Matrix H, no.uib.cipr.matrix.Vector origin,
                  no.uib.cipr.matrix.Vector dimensions) 
 {
  
   no.uib.cipr.matrix.Vector pt1 = new DenseVector(3);
  no.uib.cipr.matrix.Vector pt2 = new DenseVector(3);
   no.uib.cipr.matrix.Vector pt3 = new DenseVector(3);
  no.uib.cipr.matrix.Vector pt4 = new DenseVector(3);
  
  no.uib.cipr.matrix.Vector _pt1 = new DenseVector(3);
  no.uib.cipr.matrix.Vector _pt2 = new DenseVector(3);
   no.uib.cipr.matrix.Vector _pt3 = new DenseVector(3);
  no.uib.cipr.matrix.Vector _pt4 = new DenseVector(3);
  
   pt1.set(0,0);
   pt1.set(1,0);
   pt1.set(2,1);
   
   pt2.set(0,image.width-1);
   pt2.set(1,0);
   pt2.set(2,1);
   
   pt3.set(0,image.width-1);
   pt3.set(1,image.height-1);
   pt3.set(2,1);
   
   pt4.set(0,0);
   pt4.set(1,image.height-1);
   pt4.set(2,1);

H.mult(pt1,_pt1);
H.mult(pt2,_pt2);
H.mult(pt3,_pt3);
H.mult(pt4,_pt4);

_pt1.set(0,_pt1.get(0)/_pt1.get(2));
_pt1.set(1,_pt1.get(1)/_pt1.get(2));

_pt2.set(0,_pt1.get(0)/_pt2.get(2));
_pt2.set(1,_pt1.get(1)/_pt2.get(2));

_pt3.set(0,_pt1.get(0)/_pt3.get(2));
_pt3.set(1,_pt1.get(1)/_pt3.get(2));

_pt4.set(0,_pt1.get(0)/_pt4.get(2));
_pt4.set(1,_pt1.get(1)/_pt4.get(2));

origin.set(0,min((float)_pt1.get(0),(float)_pt4.get(0)));
origin.set(1,min((float)_pt1.get(1),(float)_pt2.get(1)));


dimensions.set(0,max((float)_pt3.get(0),(float)_pt4.get(0)));
dimensions.set(1,max((float)_pt3.get(1),(float)_pt4.get(1)));

dimensions.set(0,dimensions.get(0)-origin.get(0));
dimensions.set(1,dimensions.get(1)-origin.get(1));


}



no.uib.cipr.matrix.Vector applyHomography(Matrix homography, no.uib.cipr.matrix.Vector p){
  no.uib.cipr.matrix.Vector q = new DenseVector(new double[] {p.get(X), p.get(Y), 1.0}, false);
  no.uib.cipr.matrix.Vector r = new DenseVector(3);
  homography.mult(q, r);
  no.uib.cipr.matrix.Vector s = new DenseVector(2);
  s.set(X, r.get(X)/r.get(Z));
  s.set(Y, r.get(Y)/r.get(Z));
  return s;
}




///
/// Given a point 'p' and an image 'source' returns the bilinearly interpolated color
///
/// @param source Image from which to draw pixel values
/// @param p Point (in real coordinates) which will be interpolated
/// @pre p has (at least) 2 values
/// @returns a color value bilinearly interpolated - will be BLACK if p is not within source size
///

color getColorBilinear(PImage source, no.uib.cipr.matrix.Vector p) {
  int R, G, B;

  //  Compute the R,G,B pixel values for image coordinates p in source image
  //  using bilinear interpolation
  double x = p.get(X);
  double y = p.get(Y);
  int i = (int) Math.floor(x);
  int j = (int) Math.floor(y);
  double a = x - i;
  double b = y - j;
  
  no.uib.cipr.matrix.Vector weights;
  weights = new DenseVector(new double[] {
    (1-a)*(1-b),
       a *(1-b),
       a *   b,
    (1-a)*   b,
  });
  
  color ij = source.get(i, j);
  color iplus1j = source.get(i+1, j);
  color ijplus1 = source.get(i, j+1);
  color iplus1jplus1 = source.get(i+1, j+1);
  Matrix colors;
  colors = new DenseMatrix(new double[][] {
    { red(ij),   red(iplus1j),   red(iplus1jplus1),   red(ijplus1)   },
    { green(ij), green(iplus1j), green(iplus1jplus1), green(ijplus1) },
    { blue(ij),  blue(iplus1j),  blue(iplus1jplus1),  blue(ijplus1)  },
  });
  
  //println(colors);
  
  no.uib.cipr.matrix.Vector newc = new DenseVector(3);
  colors.mult(weights, newc);
  
  //sprintln(newColor);
  
  R = round((float) newc.get(0));
  G = round((float) newc.get(1));
  B = round((float) newc.get(2));

  return color(R,G,B);
}


PImage composite(String sourceImagePath, String sourcePointsPath, String targetImagePath,
                 String targetPointsPath) {
  // Read images and point files

  PImage sourceImage = loadImage(sourceImagePath);
  PImage targetImage = loadImage(targetImagePath);
  Matrix sourcePoints = readPoints(sourcePointsPath);
  Matrix targetPoints = readPoints(targetPointsPath);
  
  PImage maskImage;
  //if (maskImagePath != null) {
    //maskImage = loadImage(maskImagePath);
  //} else {
  maskImage = makeMask(targetImage, sourceImage, targetPoints);
  //}


  // Compute (inverse) homography matrix

  Matrix Homeo    = computeH(sourcePoints, targetPoints);
  Matrix Homeoinv = computeH(targetPoints, sourcePoints);
  
  println("Homography matrix");
  println(Homeo);
  println("Inverse homography matrix");
  println(Homeoinv);

  // Composite the images and return the resulting image
  // For this, write a function that actually does the composite and name it compositeImages.
  // This function should take the two input images (source and target image) as well as the
  // inverse homography as the input and return a composite image. You will also need to plug in
  // the the target points or a mask to determine composite region.

  PImage rectified = compositeImages(sourceImage, targetImage, maskImage, Homeoinv);
  return rectified;
}

/**
 * Create a target mask for image composition.  Returns an image the size of targetImage.
 * White within the rectangle defined by targetPoints.  Black elsewhere.
 */
PImage makeMask(PImage targetImage, PImage sourceImage, Matrix targetPoints) {
  // Assumes 4 points, specified in counter-clockwise order from top-left
  
  // Create a 100x100 white source
  int sh = sourceImage.height;
  int sw = sourceImage.width;
  PImage source = createImage(sh, sw, RGB);
  source.loadPixels();
  for (int i = 0; i < source.pixels.length; i++) {
    source.pixels[i] = color(255);
  }
  source.updatePixels();

 // Set the source points at the four corners of the source
  Matrix sourcePoints;
  sourcePoints = new DenseMatrix(new double[][] {{0, 0},{0, sh},{sw, sh},{sw, 0},});
  
  Matrix Hinv = computeH(targetPoints, sourcePoints);
  
  int h = targetImage.height;
  int w = targetImage.width;
  PImage maskImage = createImage(w, h, RGB);
  
  // Warp the source onto the mask image
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      no.uib.cipr.matrix.Vector p = new DenseVector(new double[] {x, y}, false);
      no.uib.cipr.matrix.Vector q = applyHomography(Hinv, p);
      color c = getColorBilinear(source, q);
      maskImage.set(x, y, c);
    }
  }
  
  return maskImage;
}



/**
 * Warp source onto target using Hinv homography and mask for linear blending.
 */
PImage compositeImages(PImage sourceImage, PImage targetImage, PImage maskImage, Matrix Hinv) {
  int h = targetImage.height;
  int w = targetImage.width;
  PImage result = createImage(w, h, RGB);
  
  // Warp the source onto the mask image
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      no.uib.cipr.matrix.Vector p = new DenseVector(new double[] {x, y}, false);
      no.uib.cipr.matrix.Vector q = applyHomography(Hinv, p);
      color sourceColor = getColorBilinear(sourceImage, q);
      no.uib.cipr.matrix.Vector sourceColorVec = color2Vector(sourceColor);
      no.uib.cipr.matrix.Vector targetColorVec = color2Vector(targetImage.get(x, y));
      double alpha = blue(maskImage.get(x, y))/255.0;
      // linear blending
      no.uib.cipr.matrix.Vector finalColorVec = sourceColorVec.scale(alpha).add(targetColorVec.scale(1.0-alpha));
      color c = vector2color(finalColorVec);
      result.set(x, y, c);
    }
  }
  
  return result;
}




no.uib.cipr.matrix.Vector color2Vector(color c) {
  return new DenseVector(new double[] {red(c), green(c), blue(c)}, false);
}

color vector2color(no.uib.cipr.matrix.Vector v) {
  int r = (int) v.get(0);
  int g = (int) v.get(1);
  int b = (int) v.get(2);
  return color(r, g, b);
}

//
// -- Some utility functions ----|
//


  
// general functions

Matrix readPoints(String filename) {
  String[] lines = loadStrings(filename+".txt");
  println("file is:");
  println(filename);
  DenseMatrix matrix = new DenseMatrix(lines.length, 2);
  for(int l = 0; l < lines.length; l++) {
    String[] values = splitTokens(lines[l], ", ");
    matrix.set(l,X,float(values[X]));
    matrix.set(l,Y,float(values[Y]));
  }
  return matrix;
}
//System.out.println(matrix);

Matrix computeH(Matrix src, Matrix tgt)
{
int i,j;   
int rowno;
DenseMatrix A = new DenseMatrix(2*src.numRows(),9);
// A Matrix formation
for (i=0; i<2*src.numRows(); i+=2)
{
  rowno=i/2;
  A.set(i,0, src.get(rowno, 0));
  A.set(i,0, src.get(rowno,0));        
  A.set(i,1, src.get(rowno,1));   
  A.set(i,2, 1);
  A.set(i,3, 0);
  A.set(i,4, 0);  
  A.set(i,5, 0); 
  A.set(i,6, -(src.get(rowno,0)*(tgt.get(rowno,0))));  
  A.set(i,7, -(src.get(rowno,1)*(tgt.get(rowno,0))));
  A.set(i,8, -(tgt.get(rowno,0))); 

  A.set(i+1,0, 0);
  A.set(i+1,1, 0);  
  A.set(i+1,2, 0);
  A.set(i+1,3, src.get(rowno,0));        
  A.set(i+1,4, src.get(rowno,1));   
  A.set(i+1,5, 1);
  A.set(i+1,6, -(src.get(rowno,0)*(tgt.get(rowno,1))));  
  A.set(i+1,7, -(src.get(rowno,1)*(tgt.get(rowno,1))));
  A.set(i+1,8, -(tgt.get(rowno,1)));

  
}
DenseMatrix C= new DenseMatrix(9,9);
EVD evd = new EVD(9);
double alpha=1;
A.transAmultAdd(alpha,(Matrix)A,(Matrix)C);
     try{
   evd.factor(C);
     }catch(NotConvergedException e){
    println("Error");
  }
    double eigvalues[] = evd.getRealEigenvalues();
 // println(values);
  double vect[]= new double[9];
  double w=eigvalues[0];
  int r=0;
  for(i=0;i<9;i++)
  {
    if(eigvalues[i]<=w)
    {
      w=eigvalues[i];
      r=i;
    }
  }
  
  Matrix v = evd.getLeftEigenvectors();
  for(i=0;i<9;i++)
  {
    for(j=0;j<9;j++)
    {
      if(j==r)
      {
      vect[i] = v.get(i,j);
      }
    }
  }
 
  DenseMatrix homograph = new DenseMatrix(3,3);
  for(i=0;i<3;i++)
  {
    for(j=0;j<3;j++)
    {
    homograph.set(i,j,vect[(i*3)+j]); 
    }
  }  
return homograph;
}


/// Makes small circles on the pallet image to show where the points in 'points' are
///
/// @param points points to make circles at
/// @param x horizontal offset (used for centering images, see showAndPause)
/// @param y vertical offset (used for centering images, see showAndPause)
///
void showPoints(Matrix points){
  stroke(0,0,0);
  fill(255,0,0);
  for(int r = 0; r < points.numRows(); r++) {
    ellipse((float) points.get(r,X), (float)points.get(r,Y), DOT_SIZE, DOT_SIZE);
  }
}
