package cifo;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Polygon;
import java.awt.image.BufferedImage;
import java.awt.image.PixelGrabber;
import java.util.Random;

public class Solution {

	public static final int VALUES_PER_TRIANGLE = 10;

	protected ProblemInstance instance;
	protected int[] values;
	protected double [] normalisedValues;
	protected int[] inversionOrder;
	protected int[] inversionArrayPointer;
	protected double fitness;
	protected Random r;

	public Solution(ProblemInstance instance) {
		this.instance = instance;
		r = new Random();
		initialize();
	}

	public void initialize() {
		values = new int[instance.getNumberOfTriangles() * VALUES_PER_TRIANGLE];
		normalisedValues = new double [values.length];
		inversionOrder = new int [values.length];
		inversionArrayPointer = new int [values.length];
		
		switch (Main.INITIALIZATION_METHOD) {
		case "standard": initializeStandard(); break;
		case "diverse": initializeAsDiverseAsPossible(); break;
		case "diverse_spread": initializeAsDiverseAndSpreadAsPossible(); break;
		case "big": initializeAsBigAsPossible(); break;
		default: initializeStandard(); break;
		}
	}
	


	public void initializeStandard() {
		for (int triangleIndex = 0; triangleIndex < instance.getNumberOfTriangles(); triangleIndex++) {
			// initialize HSB and Alpha
			for (int i = 0; i < 4; i++) {
				values[triangleIndex * VALUES_PER_TRIANGLE + i] = r.nextInt(256);
			}
			// initialize vertices
			for (int i = 4; i <= 8; i += 2) {
				values[triangleIndex * VALUES_PER_TRIANGLE + i] = r.nextInt(instance.getImageWidth() + 1);
				values[triangleIndex * VALUES_PER_TRIANGLE + i + 1] = r.nextInt(instance.getImageHeight() + 1);
			}
		}
	}
	
	//initializes solution with triangles as big as possible meaning, covering the whole screen. 
	//There are only 4 possible triangle initialisations possible 
	public void initializeAsBigAsPossible () {
		int xmax = instance.getImageWidth();
		int xmin = 0;
		int ymax = instance.getImageHeight();
		int ymin = 0;
		
		int [][] points = {{xmin, ymin}, {xmin, ymax}, {xmax, ymin}, {xmax, ymax}};
		
		int [][][] triangles = {{points[0], points [1], points [2]}, {points [0], points[1], points[3]}, {points[0], points [2], points [3]},{points[1], points [2], points[3]}};

		for (int triangleIndex = 0; triangleIndex < instance.getNumberOfTriangles(); triangleIndex++) {
			// initialize HSB and Alpha
			for (int j = 0; j < 4; j++) {
				values[triangleIndex*VALUES_PER_TRIANGLE + j] = r.nextInt(256);
				// System.out.println(values[triangleIndex*VALUES_PER_TRIANGLE + j]);
			}
			int k = r.nextInt(3);
			values[triangleIndex*VALUES_PER_TRIANGLE + 4] = triangles[k][0][0];
			values[triangleIndex*VALUES_PER_TRIANGLE + 5] = triangles[k][0][1];
			values[triangleIndex*VALUES_PER_TRIANGLE + 6] = triangles[k][1][0];
			values[triangleIndex*VALUES_PER_TRIANGLE + 7] = triangles[k][1][1];
			values[triangleIndex*VALUES_PER_TRIANGLE + 8] = triangles[k][2][0];
			values[triangleIndex*VALUES_PER_TRIANGLE + 9] = triangles[k][2][1];
		}
	}

	private void initializeAsDiverseAndSpreadAsPossible() {
		//100 Triangles * 3 points per triangle  = 300 points
		//200*200 field
		int [][] colorValues = initializeDiverseColorValues();
		int totalNumberPoints = instance.getNumberOfTriangles()*3;
		double pointPerAxis = Math.sqrt(totalNumberPoints);
		int pointPerXAxis = (int)Math.ceil(pointPerAxis);
		int pointPerYAxis = (int)Math.floor(pointPerAxis);
		double numbersPerXAxis = instance.getImageWidth()/(pointPerXAxis-1.0);
		double numbersPerYAxis = instance.getImageHeight()/(pointPerYAxis-1.0);
		int [][][] verticeValues = new int [pointPerXAxis][pointPerYAxis][2];
		int currentNumberPoints = 0;
		for (int i=0; i<pointPerXAxis; i++) {
			int xVertice = (int)(i*numbersPerXAxis); 
			for (int j = 0; j<pointPerYAxis; j++) {
				if (currentNumberPoints<totalNumberPoints) {
					verticeValues [i][j][0] = xVertice;
					verticeValues [i][j][1] = (int)(j*numbersPerYAxis);
					currentNumberPoints++;
				} else {
					verticeValues [i][j] = null;
				}
			}
		}
		int randomPointX;
		int randomPointY;
		boolean solutionIsNew;
		for (int triangleIndex = 0; triangleIndex < instance.getNumberOfTriangles(); triangleIndex++) {
			// initialize HSB and Alpha
			for (int j = 0; j < 4; j++) {
				solutionIsNew=false;
				randomPointX = r.nextInt(instance.getNumberOfTriangles());
				while (!solutionIsNew) {
					if (colorValues[j][randomPointX]!=-1) {
						values[triangleIndex*VALUES_PER_TRIANGLE + j] = colorValues[j][randomPointX];
						colorValues[j][randomPointX] = -1;
						solutionIsNew=true;
					}
					randomPointX = r.nextInt(instance.getNumberOfTriangles());
				}
			}
			randomPointX = r.nextInt(pointPerXAxis);
			randomPointY = r.nextInt(pointPerYAxis);
			for (int j=0; j<5; j = j+2) {
				solutionIsNew=false;
				randomPointX = r.nextInt(pointPerXAxis);
				randomPointY = r.nextInt(pointPerYAxis);

				while (!solutionIsNew) {
					if (verticeValues[randomPointX][randomPointY]!=null) {
						values[triangleIndex*VALUES_PER_TRIANGLE + j + 4] = verticeValues[randomPointX][randomPointY][j%2];
						values[triangleIndex*VALUES_PER_TRIANGLE + j+1 + 4] = verticeValues[randomPointX][randomPointY][(j+1)%2];
						verticeValues[randomPointX][randomPointY] = null;
						solutionIsNew=true;
					}
					randomPointX = r.nextInt(pointPerXAxis);
					randomPointY = r.nextInt(pointPerYAxis);
				}
			}
		}
		
	}
	
	public int [] [] initializeDiverseColorValues() {
		int [][] colorValues = new int [4][instance.getNumberOfTriangles()];
		for (int i=0; i<instance.getNumberOfTriangles(); i++) {
			colorValues[0][i] = 255 - i*255/(instance.getNumberOfTriangles());
			colorValues[1][i] = colorValues[0][i] - (int)Math.round(1.0/3 * ((i+1)*255.0/instance.getNumberOfTriangles()-(i*255.0/instance.getNumberOfTriangles())));
			colorValues[2][i] = colorValues[0][i] - (int)Math.round(2.0/3 * ((i+1)*255.0/instance.getNumberOfTriangles()-(i*255.0/instance.getNumberOfTriangles())));
			colorValues[3][i] = colorValues[0][i] - (int)Math.round(3.0/3 * ((i+1)*255.0/instance.getNumberOfTriangles()-(i*255.0/instance.getNumberOfTriangles())));
			}
		return colorValues;
	}
	
	public int [] [] initializeDiverseVerticeValues() {
		int [][] verticeValues = new int [6][instance.getNumberOfTriangles()];
		for (int i=0; i<instance.getNumberOfTriangles(); i++) {
			verticeValues[0][i] = verticeValues[2][i] = verticeValues[4][i] = instance.getImageWidth()-1 - i*(instance.getImageWidth()-1)/(instance.getNumberOfTriangles()-1);
			verticeValues[1][i] = verticeValues[3][i] = verticeValues[5][i] = instance.getImageHeight()-1 - i*(instance.getImageHeight()-1)/(instance.getNumberOfTriangles()-1);
		}
		return verticeValues;
	}
	
	public void initializeAsDiverseAsPossible () {			
		int [][] colorValues;
		int [][] verticeValues;
		
		colorValues = this.initializeDiverseColorValues();
		verticeValues = this.initializeDiverseVerticeValues();
		int randomPoint = r.nextInt(instance.getNumberOfTriangles());
		boolean solutionIsNew;
		for (int triangleIndex = 0; triangleIndex < instance.getNumberOfTriangles(); triangleIndex++) {
			// initialize HSB and Alpha
			for (int j = 0; j < 4; j++) {
				solutionIsNew=false;
				while (!solutionIsNew) {
					if (colorValues[j][randomPoint]!=-1) {
						values[triangleIndex*VALUES_PER_TRIANGLE + j] = colorValues[j][randomPoint];
						colorValues[j][randomPoint] = -1;
						solutionIsNew=true;
					}
					randomPoint = r.nextInt(instance.getNumberOfTriangles());
				}
			}
			for (int j=0; j<6; j++) {
				solutionIsNew=false;
				while (!solutionIsNew) {
					if (verticeValues[j][randomPoint]!=-1) {
						values[triangleIndex*VALUES_PER_TRIANGLE + j+4] = verticeValues[j][randomPoint];
						verticeValues[j][randomPoint] = -1;
						solutionIsNew=true;
					}
					randomPoint = r.nextInt(instance.getNumberOfTriangles());
				}
			}
		}

	}	
	
	//initialize with pre-given triangle values
	public void initialize (int [] values) {
		this.values = values;
	}
	
	//values [] has to be of size triangleNumber*4
	public void initializeColorValues(int [] values) {
		this.values = new int[instance.getNumberOfTriangles() * VALUES_PER_TRIANGLE];
		
		for (int triangleIndex = 0; triangleIndex < instance.getNumberOfTriangles(); triangleIndex++) {
			for (int i = 0; i < 4; i++) {
				this.values[triangleIndex * VALUES_PER_TRIANGLE + i] = values [i];
				this.values[triangleIndex * VALUES_PER_TRIANGLE + i] = values [triangleIndex*4 + i];
			}
		}
	}
	
	//values [] has to be of size triangleNumber*6
	public void initializeVerticeValues(int [] values) {
		this.values = new int[instance.getNumberOfTriangles() * VALUES_PER_TRIANGLE];
		
		for (int triangleIndex = 0; triangleIndex < instance.getNumberOfTriangles(); triangleIndex++) {
			for (int i = 4; i < 10; i++) {
				this.values[triangleIndex * VALUES_PER_TRIANGLE + i] = values [i-4];
				this.values[triangleIndex * VALUES_PER_TRIANGLE + i] = values [triangleIndex * 6 + i-4];
			}
		}
	}

	public void evaluate() {
		BufferedImage generatedImage = createImage();
		int[] generatedPixels = new int[generatedImage.getWidth() * generatedImage.getHeight()];
		PixelGrabber pg = new PixelGrabber(generatedImage, 0, 0, generatedImage.getWidth(), generatedImage.getHeight(),
				generatedPixels, 0, generatedImage.getWidth());
		try {
			pg.grabPixels();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		int[] targetPixels = instance.getTargetPixels();
		long sum = 0;
		for (int i = 0; i < targetPixels.length; i++) {
			int c1 = targetPixels[i];
			int c2 = generatedPixels[i];
			int red = ((c1 >> 16) & 0xff) - ((c2 >> 16) & 0xff);
			int green = ((c1 >> 8) & 0xff) - ((c2 >> 8) & 0xff);
			int blue = (c1 & 0xff) - (c2 & 0xff);
			sum += red * red + green * green + blue * blue;
		}

		this.fitness = Math.sqrt(sum);
	}

	public Solution applyMutation(int triangleIndex) {
		Solution temp = this.copy();
		int valueIndex = r.nextInt(VALUES_PER_TRIANGLE);
		if (valueIndex < 4) {
			temp.values[triangleIndex * VALUES_PER_TRIANGLE + valueIndex] = r.nextInt(256);
		} else {
			if (valueIndex % 2 == 0) {
				temp.values[triangleIndex * VALUES_PER_TRIANGLE + valueIndex] = r.nextInt(instance.getImageWidth() + 1);
			} else {
				temp.values[triangleIndex * VALUES_PER_TRIANGLE + valueIndex] = r
						.nextInt(instance.getImageHeight() + 1);
			}
		}
		temp.printValues();
		return temp;
	}

	public void draw() {
		BufferedImage generatedImage = createImage();
		Graphics g = ProblemInstance.view.getFittestDrawingView().getMainPanel().getGraphics();
		g.drawImage(generatedImage, 0, 0, ProblemInstance.view.getFittestDrawingView());
	}

	public void print() {
		System.out.printf("Fitness: %.1f\n", fitness);
	}

	public int getValue(int index) {
		return values[index];
	}

	public void setValue(int index, int value) {
		values[index] = value;
	}

	public int getHue(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 0];
	}

	public int getSaturation(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 1];
	}

	public int getBrightness(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 2];
	}

	public int getAlpha(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 3];
	}

	public int getXFromVertex1(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 4];
	}

	public int getYFromVertex1(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 5];
	}

	public int getXFromVertex2(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 6];
	}

	public int getYFromVertex2(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 7];
	}

	public int getXFromVertex3(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 8];
	}

	public int getYFromVertex3(int triangleIndex) {
		return values[triangleIndex * VALUES_PER_TRIANGLE + 9];
	}

	public void setHue(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 0] = value;
	}

	public void setSaturation(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 1] = value;
	}

	public void setBrightness(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 2] = value;
	}

	public void setAlpha(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 3] = value;
	}

	public void setXFromVertex1(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 4] = value;
	}

	public void setYFromVertex1(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 5] = value;
	}

	public void setXFromVertex2(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 6] = value;
	}

	public void setYFromVertex2(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 7] = value;
	}

	public void setXFromVertex3(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 8] = value;
	}

	public void setYFromVertex3(int triangleIndex, int value) {
		values[triangleIndex * VALUES_PER_TRIANGLE + 9] = value;
	}

	public int[] getVertex1(int triangleIndex) {
		return new int[] { getXFromVertex1(triangleIndex), getYFromVertex1(triangleIndex) };
	}

	public int[] getVertex2(int triangleIndex) {
		return new int[] { getXFromVertex2(triangleIndex), getYFromVertex2(triangleIndex) };
	}

	public int[] getVertex3(int triangleIndex) {
		return new int[] { getXFromVertex3(triangleIndex), getYFromVertex3(triangleIndex) };
	}

	public ProblemInstance getInstance() {
		return instance;
	}

	public int[] getValues() {
		return values;
	}

	public double getFitness() {
		return fitness;
	}
	
	public String printValues () {
		String print = "";
		for (int i=0; i<getValues().length-1; i++) {
			print += getValues()[i] + ";";
		}
		return print += getValues()[getValues().length-1];
	}

	public Solution copy() {
		Solution temp = new Solution(instance);
		for (int i = 0; i < values.length; i++) {
			temp.values[i] = values[i];
		}
		temp.fitness = fitness;
		return temp;
	}

	private BufferedImage createImage() {
		BufferedImage target = instance.getTargetImage();
		BufferedImage generatedImage = new BufferedImage(target.getWidth(), target.getHeight(),
				BufferedImage.TYPE_INT_ARGB);
		Graphics generatedGraphics = generatedImage.getGraphics();

		generatedGraphics.setColor(Color.GRAY);
		generatedGraphics.fillRect(0, 0, generatedImage.getWidth(), generatedImage.getHeight());
		for (int triangleIndex = 0; triangleIndex < instance.getNumberOfTriangles(); triangleIndex++) {
			generatedGraphics.setColor(expressColor(triangleIndex));
			generatedGraphics.fillPolygon(expressPolygon(triangleIndex));
		}
		return generatedImage;
	}

	private Color expressColor(int triangleIndex) {
		int hue = getHue(triangleIndex);
		int saturation = getSaturation(triangleIndex);
		int brightness = getBrightness(triangleIndex);
		int alpha = getAlpha(triangleIndex);
		Color c = Color.getHSBColor(hue / 255.0f, saturation / 255.0f, brightness / 255.0f);
		return new Color(c.getRed(), c.getGreen(), c.getBlue(), alpha);
	}

	private Polygon expressPolygon(int triangleIndex) {
		int[] xs = new int[] { getXFromVertex1(triangleIndex), getXFromVertex2(triangleIndex),
				getXFromVertex3(triangleIndex) };
		int[] ys = new int[] { getYFromVertex1(triangleIndex), getYFromVertex2(triangleIndex),
				getYFromVertex3(triangleIndex) };
		return new Polygon(xs, ys, 3);
	}

	public void applyInversion() {	
		int [] newValues = values.clone();

		int randomValue = 0;
		for (int i=0; i<instance.getNumberOfTriangles()*VALUES_PER_TRIANGLE; i++) {
			randomValue = r.nextInt(instance.getNumberOfTriangles()*VALUES_PER_TRIANGLE);
			boolean applicable = false;
			if ((randomValue)%10<4 && i%10<4) applicable = true;
			if ((randomValue)%10>=4 && i%10>=4) applicable = true;
			
			while (!isEmpty(randomValue) || !applicable) {
				randomValue = r.nextInt(instance.getNumberOfTriangles()*VALUES_PER_TRIANGLE) ;
				if ((randomValue)%10<4 && i%10<4) applicable = true;
				if ((randomValue)%10>=4 && i%10>=4) applicable = true;
			}
			values[randomValue] = newValues[i];
			inversionOrder[randomValue] = i+1;
			inversionArrayPointer[i] = randomValue;
		}
		for (int i=0; i<inversionOrder.length; i++) {
			inversionOrder[i] -= 1;
		}
	}
	
	public void applyNormalisedInversion() {
		double [] newValues = normalise(values);
		//reinitialisation to make sure they are empty
		
		//int [] newValues = values.clone();
		
//		inversionOrder = new int [values.length];
//		inversionArrayPointer = new int [values.length];
		int randomValue = 0;
		for (int i=0; i<instance.getNumberOfTriangles()*VALUES_PER_TRIANGLE; i++) {
			randomValue = r.nextInt(instance.getNumberOfTriangles()*VALUES_PER_TRIANGLE) + 1;
			
			while (!isEmpty(randomValue)) {
				randomValue = r.nextInt(instance.getNumberOfTriangles()*VALUES_PER_TRIANGLE) + 1;
			}
			normalisedValues[randomValue-1] = newValues[i];
			inversionOrder[i] = randomValue;
			inversionArrayPointer[randomValue-1] = i;
		}
		for (int i=0; i<inversionOrder.length; i++) {
			inversionOrder[i] -= 1;
		}
	}
	
	public void applyReordering() {
		int [] newValues = values.clone();
		for (int i=0; i<values.length; i++) {
			values[inversionOrder[i]] = newValues[i];
		}
	}
	
	public void applyNormalisedReordering() {
		double [] newValues = normalisedValues.clone();
		for (int i=0; i<values.length; i++) {
			normalisedValues[i] = newValues[inversionOrder[i]];
		}
		values = denormalise(normalisedValues);
	}

	private double[] normalise(int [] startValues) {
		double [] normalisedValues = new double [startValues.length];
		for (int i=0; i<instance.getNumberOfTriangles(); i++) {
			for (int j=0; j<4; j++) {
				normalisedValues [i*VALUES_PER_TRIANGLE + j] = (startValues[i*VALUES_PER_TRIANGLE + j]/256.0) * 1000000;
			}
			for (int j=4; j<9; j=j+2) {
				normalisedValues [i*VALUES_PER_TRIANGLE + j] = (startValues [i*VALUES_PER_TRIANGLE + j]/(double)instance.getImageWidth()) * 1000000;
				normalisedValues [i*VALUES_PER_TRIANGLE + j+1] = (startValues [i*VALUES_PER_TRIANGLE + j+1]/(double)instance.getImageHeight() * 1000000);
			}
		}
		return normalisedValues;
	}
	
private int[] denormalise(double[] normalisedValues) {
//		for (int j=0; j<instance.getNumberOfTriangles()*VALUES_PER_TRIANGLE; j++) {
//			System.out.println((startValues [j]));
//		}
	
	int [] denormalisedValues = new int [normalisedValues.length];
		
		for (int i=0; i<instance.getNumberOfTriangles(); i++) {
			for (int j=0; j<4; j++) {
				denormalisedValues [i*VALUES_PER_TRIANGLE + j] = (int)Math.floor((normalisedValues[i*VALUES_PER_TRIANGLE + j]/ 1000000.0) * 255);
			}
			for (int j=4; j<9; j=j+2) {
				denormalisedValues [i*VALUES_PER_TRIANGLE + j] = (int)Math.floor((normalisedValues[i*VALUES_PER_TRIANGLE + j] / 1000000.0) * instance.getImageWidth());
				denormalisedValues [i*VALUES_PER_TRIANGLE + j+1] = (int)Math.floor((normalisedValues[i*VALUES_PER_TRIANGLE + j+1] / 1000000.0) * instance.getImageHeight());
			}
		}
//		for (int j=0; j<instance.getNumberOfTriangles()*VALUES_PER_TRIANGLE; j++) {
//			if (startValues [j]>256) {
//				System.out.println((j) + " : " + (startValues [j]));
//			}
//		}
		return denormalisedValues;
	}

	private boolean isEmpty(int randomValue) {
		for (int i=0; i<inversionOrder.length; i++) {
			if (inversionOrder[randomValue]>0) return false;
		}
		return true;
	}
	
	public int[] getInversionOrder() {
		return inversionOrder;
	}

	public int[] getInversionArrayPointer() {
		return inversionArrayPointer;
	}
	
	public int getInversionOrderIndex(int index) {
		return inversionOrder[index];
	}
	
	public void setInversionOrderIndex(int index, int value) {
		this.inversionOrder[index] = value;
	}
	
	public int getInversionArrayPointerIndex(int index) {
		return inversionArrayPointer[index];
	}
	
	public double[] getNormalisedValues() {
		return normalisedValues;
	}
	
	public double getNormalisedValue (int index) {
		return this.normalisedValues[index];
	}
	
	public void setNormalisedValue (int index, double value) {
		this.normalisedValues[index] = value;
	}
}