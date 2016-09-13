package SVTFactory;

import java.io.IOException;
import java.io.Writer;

import org.jlab.geom.base.ConstantProvider;
import org.jlab.geom.geant.Geant4Basic;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Transformation3D;
import org.jlab.geom.prim.Triangle3D;

import Alignment.AlignmentFactory;
import Misc.Util;

/**
 * <h1> Geometry for the SVT </h1>
 * 
 * Processes fiducial survey data into alignment shifts, and applies those shifts to a given point or volume.
 * 
 * @author pdavies
 * @version 0.2.2
 */
public class SVTAlignmentFactory
{
	private static String filenameSurveyIdeal;
	private static String filenameSurveyMeasured;
	
	//private static String filenameDistances = "measured_distances.dat";
	//private static Writer outputDistances;
	
	//private static String filenameIdealFiducials = "survey_Ideal_reformat2.dat";
	//private static Writer outputIdealFiducials;
	
	//private static String filenameMeasuredFiducials = "survey_measured_reformat2.dat";
	//private static Writer outputMeasuredFiducials;
	
	private static double[][] dataSurveyIdeal, dataSurveyMeasured;
	
	/**
	 * Constructs a new geometry factory for sensor strips.
	 * Please run {@code SVTConstants.connect() } first.
	 * 
	 * @param cp a DatabaseConstantProvider that has loaded the necessary tables
	 * @param aInputSurveyIdeal a filename for the ideal data of the fiducial survey
	 * @param aInputSurveyMeasured a filename for the measured data of the fiducial survey
	 */
	public static void setup( ConstantProvider cp, String aInputSurveyIdeal, String aInputSurveyMeasured )
	{
		SVTConstants.load( cp );
		filenameSurveyIdeal = aInputSurveyIdeal;
		filenameSurveyMeasured = aInputSurveyMeasured;
		
		try
		{
			dataSurveyIdeal = Util.inputTaggedData( filenameSurveyIdeal, 3 ); // RSF (X Y Z)
			dataSurveyMeasured = Util.inputTaggedData( filenameSurveyMeasured, 3 ); // RSF (X Y Z)
		}
		catch( IOException e ){ e.printStackTrace(); }
		
		if( dataSurveyIdeal == null || dataSurveyMeasured == null )
			throw new IllegalArgumentException("no data");
	}

	
	
	/**
	 * Calculates the alignment shifts between two sets of fiducial data, and writes them to the given file.
	 * 
	 * @param aDataIdeal fiducial data
	 * @param aDataMeasured fiducial data
	 * @param aOutputFile a filename
	 * @return double[][] translations and axis-angle rotations of the form { tx, ty, tz, rx, ry, rz, ra } relative to the first data set
	 */
	public static double[][] calcShifts( double[][] aDataIdeal, double[][] aDataMeasured, String aOutputFile )
	{
		String outputLine; Writer outputShifts;
		double[][] dataShifts = AlignmentFactory.calcShifts( SVTConstants.NTOTALSECTORS, aDataIdeal, aDataMeasured );
		
		outputShifts = Util.openOutputDataFile( aOutputFile );	
		
		for( int l = 0; l < SVTConstants.NTOTALSECTORS; l++ )
		{
			int[] rs = SVTConstants.convertIndex2RegionSector( l );
			String fmt = "R%dS%02d % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f\n";
			outputLine = String.format(fmt, rs[0]+1, rs[1]+1,
					dataShifts[l][0], dataShifts[l][1], dataShifts[l][2], dataShifts[l][3], dataShifts[l][4], dataShifts[l][5], Math.toDegrees(dataShifts[l][6]) );
			//outputLine = String.format("R%dS%02d % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f\n", region+1, sector+1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
			Util.writeLine( outputShifts, outputLine );
			//System.out.print(outputLine);
		}
		Util.closeOutputDataFile( aOutputFile, outputShifts );
		return dataShifts;
	}
	
	
	/**
	 * Calculates the difference between two sets of fiducial data, and writes them to the given file. 
	 * 
	 * @param aDataIdeal first set of fiducial data
	 * @param aDataMeasured second set of fiducial data
	 * @param aOutputFile a filename
	 * @return double[][] point differences relative to the first data set
	 */
	public static double[][] calcDeltas( double[][] aDataIdeal, double[][] aDataMeasured, String aOutputFile )
	{
		String outputLine; Writer outputWriter;
		outputWriter = Util.openOutputDataFile( aOutputFile );
		
		double[][] dataDeltasMeasuredFromIdeal = AlignmentFactory.calcDeltas( SVTConstants.NTOTALSECTORS*SVTConstants.NFIDUCIALS, 3, aDataIdeal, aDataMeasured );

		for( int k = 0; k < SVTConstants.NTOTALSECTORS*SVTConstants.NFIDUCIALS; k++ )
		{
			int[] rsf = SVTConstants.convertIndex2RegionSectorFiducial( k );
			
			double[] data = dataDeltasMeasuredFromIdeal[k];
			double radiusSpherical = Math.sqrt( Math.pow(data[0], 2 ) + Math.pow(data[1], 2 ) + Math.pow(data[2], 2 ) );
			double radiusCylindrical = Math.sqrt( Math.pow(data[0], 2 ) + Math.pow(data[1], 2 ) );
			
			outputLine = String.format("R%dS%02dF%d % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f\n", rsf[0]+1, rsf[1]+1, rsf[2]+1,
					data[0], data[1], data[2], radiusCylindrical, radiusSpherical );
			Util.writeLine( outputWriter, outputLine );
			//System.out.print(outputLine);
		}
		Util.closeOutputDataFile( aOutputFile, outputWriter );
		return dataDeltasMeasuredFromIdeal;
	}
	
	
	/**
	 * Calculates the difference between two sets of distance data, and writes them to the given file.
	 * 
	 * @param aDataIdeal first set of distance data
	 * @param aDataMeasured second set of distance data
	 * @param aOutputFile a filename
	 * @return double[][][] distance differences relative to the first data set in {value,uncertainty} pairs
	 */
	public static double[][][] calcDistanceDeltas( double[][][] aDataIdeal, double[][][] aDataMeasured, String aOutputFile )
	{
		String outputLine; Writer outputWriter;
		outputWriter = Util.openOutputDataFile( aOutputFile );
		
		double[][][] data = AlignmentFactory.calcDeltas( SVTConstants.NTOTALSECTORS, 3, 2, aDataIdeal, aDataMeasured );
		
		for( int j = 0; j < SVTConstants.NTOTALSECTORS; j++ )
		{
			int[] rs = SVTConstants.convertIndex2RegionSector( j );
			outputLine = String.format("R%dS%02d % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f\n", rs[0]+1, rs[1]+1,
										data[j][0][0], data[j][1][0], data[j][2][0], 
										data[j][0][1], data[j][1][1], data[j][2][1] );
			Util.writeLine( outputWriter, outputLine );
		}
		Util.closeOutputDataFile( aOutputFile, outputWriter );
		return data;
	}
	
	
	/**
	 * Calculates the distances between points for each triangle in a set of fiducial data, and writes them to the given file.
	 * 
	 * @param aData the set of fiducial data
	 * @param aUncertainty uncertainty in the measurement of the coordinates
	 * @param aOutputFile a filename
	 * @return double[][][] distances in {value,uncertainty} pairs
	 */
	public static double[][][] calcTriangleSides( double[][] aData, double aUncertainty, String aOutputFile )
	{
		String outputLine; Writer outputWriter;
		outputWriter = Util.openOutputDataFile( aOutputFile );
		
		double[][][] distances = new double[SVTConstants.NTOTALSECTORS][][];
		
		for( int j = 0; j < SVTConstants.NTOTALSECTORS; j++ )
		{
			int[] rs = SVTConstants.convertIndex2RegionSector( j );
			int k0 = SVTConstants.convertRegionSectorFiducial2Index(rs[0], rs[1], 0 );
			int k1 = SVTConstants.convertRegionSectorFiducial2Index(rs[0], rs[1], 1 );
			int k2 = SVTConstants.convertRegionSectorFiducial2Index(rs[0], rs[1], 2 );
			//System.out.printf("j%2d r%d s%2d k%3d %3d %3d\n",j, rs[0], rs[1], k0, k1, k2 );
			
			Point3D[] pos3Ds = new Point3D[]{ Util.toPoint3D(aData[k0]), Util.toPoint3D(aData[k1]), Util.toPoint3D(aData[k2]) };
			
			distances[j] = new double[][]{ Util.calcDistance( pos3Ds[0], pos3Ds[1], aUncertainty, aUncertainty ),
										   Util.calcDistance( pos3Ds[1], pos3Ds[2], aUncertainty, aUncertainty ),
										   Util.calcDistance( pos3Ds[2], pos3Ds[0], aUncertainty, aUncertainty ) };
			
			outputLine = String.format("R%dS%02d % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f % 8.3f\n", rs[0]+1, rs[1]+1,
										distances[j][0][0], distances[j][1][0], distances[j][2][0], 
										distances[j][0][1], distances[j][1][1], distances[j][2][1] );
			
			Util.writeLine( outputWriter, outputLine );
		}
		
		Util.closeOutputDataFile( aOutputFile, outputWriter );
		return distances;
	}
	
	
	/**
	 * Applies the given alignment shift to the given point.
	 * 
	 * @param aPos a volume in the lab frame
	 * @param aShift a translation and axis-angle rotation of the form { tx, ty, tz, rx, ry, rz, ra }
	 * @param aNominalCenter a point about which to rotate the first point (for example the midpoint of the nominal fiducials)
	 */
	public static void applyShift( Point3D aPos, double[] aShift, Point3D aNominalCenter )
	{
		AlignmentFactory.applyShift( aPos, aShift, aNominalCenter, 1.0, 1.0 );
	}
	
	
	/**
	 * Applies the given alignment shift to the given volume.
	 * 
	 * @param aVol a volume in the lab frame
	 * @param aShift a translation and axis-angle rotation of the form { tx, ty, tz, rx, ry, rz, ra }
	 * @param aNominalCenter a point about which to rotate the first point (for example the midpoint of the nominal fiducials)
	 */
	public static void applyShift( Geant4Basic aVol, double[] aShift, Point3D aNominalCenter )
	{
		AlignmentFactory.applyShift( aVol, aShift, aNominalCenter, 1.0, 1.0 );
	}
	
	
	/**
	 * Returns locations of shifted fiducial points.
	 * 
	 * @return double[][] an array of data in fiducial survey format.
	 */
	public static double[][] getShiftedFiducialData()
	{
		double [][] data = new double[SVTConstants.NTOTALSECTORS*SVTConstants.NFIDUCIALS][];
		for( int region = 0; region < SVTConstants.NREGIONS; region++ )
			for( int sector = 0; sector < SVTConstants.NSECTORS[region]; sector++ )
			{
				Point3D fidPos3Ds[] = getShiftedFiducials( region, sector );
				for( int fid = 0; fid < SVTConstants.NFIDUCIALS; fid++ )
					data[SVTConstants.convertRegionSectorFiducial2Index( region, sector, fid )] = new double[]{ fidPos3Ds[fid].x(), fidPos3Ds[fid].y(), fidPos3Ds[fid].z() };
			}
		return data;
	}
	
	
	/**
	 * Returns locations of nominal fiducial points.
	 * 
	 * @return double[][] an array of data in fiducial survey format.
	 */
	public static double[][] getFactoryIdealFiducialData()
	{
		double [][] data = new double[SVTConstants.NTOTALSECTORS*SVTConstants.NFIDUCIALS][];
		for( int region = 0; region < SVTConstants.NREGIONS; region++ )
			for( int sector = 0; sector < SVTConstants.NSECTORS[region]; sector++ )
			{
				Point3D fidPos3Ds[] = getIdealFiducials( region, sector );
				for( int fid = 0; fid < SVTConstants.NFIDUCIALS; fid++ )
				{
					data[SVTConstants.convertRegionSectorFiducial2Index( region, sector, fid )] 
							= new double[]{ fidPos3Ds[fid].x(), fidPos3Ds[fid].y(), fidPos3Ds[fid].z() };
				}
			}
		return data;
	}
	
	
	/**
	 * Returns a set of fiducial points for a sector module after the alignment shifts been applied.
	 * These indices start from 0.
	 * 
	 * @param aRegion an index starting from 0
	 * @param aSector an index starting from 0
	 * @return Point3D[] an array of fiducial points in the order Cu+, Cu-, Pk
	 * @throws IllegalArgumentException indices out of bounds
	 */
	public static Point3D[] getShiftedFiducials( int aRegion, int aSector ) throws IllegalArgumentException
	{
		if( aRegion < 0 || aRegion > SVTConstants.NREGIONS-1 ){ throw new IllegalArgumentException("region out of bounds"); }
		if( aSector < 0 || aSector > SVTConstants.NSECTORS[aRegion]-1 ){ throw new IllegalArgumentException("sector out of bounds"); }
		
		Point3D[] fidPos3Ds = getIdealFiducials( aRegion, aSector ); // lab frame
		Triangle3D fidTri3D = new Triangle3D( fidPos3Ds[0], fidPos3Ds[1], fidPos3Ds[2] );
		
		for( int f = 0; f < SVTConstants.NFIDUCIALS; f++ )
			applyShift( fidPos3Ds[f], SVTConstants.getAlignmentShiftData()[SVTConstants.convertRegionSector2Index( aRegion, aSector )], fidTri3D.center() );
		
		return fidPos3Ds;
	}
	
	
	/**
	 * Returns a set of fiducial points for a sector module before any alignment shifts been applied.
	 * 
	 * @param aRegion an index starting from 0
	 * @param aSector an index starting from 0
	 * @return Point3D[] an array of fiducial points in the order Cu+, Cu-, Pk
	 * @throws IllegalArgumentException indices out of bounds
	 */
	public static Point3D[] getIdealFiducials( int aRegion, int aSector ) throws IllegalArgumentException // lab frame
	{
		if( aRegion < 0 || aRegion > SVTConstants.NREGIONS-1 ){ throw new IllegalArgumentException("region out of bounds"); }
		if( aSector < 0 || aSector > SVTConstants.NSECTORS[aRegion]-1 ){ throw new IllegalArgumentException("sector out of bounds"); }
		
		Point3D[] fidPos3Ds = new Point3D[] { createFiducial(0), createFiducial(1), createFiducial(2) }; // relative to fiducial origin
		
		double fidOriginZ = SVTConstants.Z0ACTIVE[aRegion] - SVTConstants.DEADZNLEN - SVTConstants.FIDORIGINZ;
		double copperWideThk = 2.880;
		double radius = SVTConstants.SUPPORTRADIUS[aRegion] + copperWideThk;
		
		Transformation3D labFrame = SVTConstants.getLabFrame( aRegion, aSector, radius, fidOriginZ );
		
		for( int f = 0; f < SVTConstants.NFIDUCIALS; f++ )
			labFrame.apply( fidPos3Ds[f] );
		
		return fidPos3Ds;
	}
	
	/**
	 * Returns the fiducial center for a sector module before any alignment shifts been applied.
	 * 
	 * @param aRegion an index starting from 0
	 * @param aSector an index starting from 0
	 * @return Point3D the mean average point of the 3 fiducial points (Cu+, Cu-, Pk)
	 * @throws IllegalArgumentException indices out of bounds
	 */
	public static Point3D getIdealFiducialCenter( int aRegion, int aSector ) throws IllegalArgumentException
	{
		Point3D[] fidPos3Ds = getIdealFiducials( aRegion, aSector );
		return new Triangle3D( fidPos3Ds[0], fidPos3Ds[1], fidPos3Ds[2] ).center();
	}
	
	/**
	 * Returns a fiducial point on a sector module in the local frame.
	 * 
	 * @param aFid an index for the desired point: 0, 1, 2
	 * @return Point3D one of 3 fiducial points: Cu+, Cu-, Pk
	 * @throws IllegalArgumentException indices out of bounds
	 */
	public static Point3D createFiducial( int aFid ) throws IllegalArgumentException // local frame
	{
		if( aFid < 0 || aFid > SVTConstants.NFIDUCIALS-1 ){ throw new IllegalArgumentException("region out of bounds"); }
		
		Point3D fidPos = new Point3D();
		
		switch( aFid )
		{
		case 0: // Cu +
			fidPos.set( SVTConstants.FIDCUX, 0.0, -SVTConstants.FIDCUZ );
			break;
		case 1: // Cu -
			fidPos.set( -SVTConstants.FIDCUX, 0.0, -SVTConstants.FIDCUZ );
			break;
		case 2: // Pk
			fidPos.set( -SVTConstants.FIDPKX, 0.0, SVTConstants.FIDPKZ0 + SVTConstants.FIDPKZ1 );
			break;
		}
		return fidPos;
	}
	
	
	
	public static double[][] getDataSurveyIdeal()
	{
		return dataSurveyIdeal;
	}
	
	
	
	public static double[][] getDataSurveyMeasured()
	{
		return dataSurveyMeasured;
	}
}

