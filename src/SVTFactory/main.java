package SVTFactory;

import java.io.Writer;

//import org.jlab.clasrec.utils.DatabaseConstantProvider; // 2.4
import org.jlab.detector.calib.utils.DatabaseConstantProvider; // 3.0
import org.jlab.geom.geant.Geant4Basic;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Triangle3D;
import org.jlab.geom.prim.Vector3D;

import Misc.*;
import VolumeExporter.*;

public class main
{
	public static void main(String[] args)
	{
		//ConstantProvider myconstants = DataBaseLoader.getGeometryConstants( DetectorType.BST );
		
		//FTOFFactory myfactory = new FTOFFactory();
		//FTOFDetector mydetector = myfactory.createDetectorCLAS(myconstants);
		//System.exit(0);
		
		//BSTFactory myfactory = new BSTFactory();
		//BSTDetector mydetector = myfactory.createDetectorCLAS(myconstants);
		//System.exit(0);
		
		/*ConstantsManager cm = new ConstantsManager();
		cm.init(Arrays.asList(new String[]{ "/geometry/cvt/svt" }));
		IndexedTable bstTable = cm.getConstants( 10, "/geometry/cvt/svt");*/
		//System.exit(0);
		
		/*Vector3D axis = new Vector3D( 0, 1, 1 ).asUnit();
		double angle = 180;
		double[] aa = new double[]{ axis.x(), axis.y(), axis.z(), Math.toRadians( angle ) };
		Matrix maa = Matrix.convertRotationAxisAngleToMatrix( aa );
		
		Matrix v0 = new Matrix( 3, 1, new double[]{ 0, 0, 1 } );
		Matrix v1 = Matrix.matMul( maa, v0 );
		
		System.out.print("axis ");
		axis.show();
		System.out.println("angle "+ angle );
		maa.show("axis angle matrix");
		v0.show("vector before");
		v1.show("vector after");*/
		
		/*Geant4Basic topVol = new Geant4Basic("top", "Box", 1, 1, 2 );
		
		Geant4Basic myVol = new Geant4Basic("myVol", "Box", 2, 2, 4 );
		myVol.setMother( topVol );
		myVol.setPosition( 0, -5, 0 );
		//myVol.setRotation("xyz", 0, -Math.toRadians(30), 0 );
		myVol.setRotation("xyz", 0, 0, -Math.toRadians(-45) );
		
		Geant4Basic yourVol = new Geant4Basic("yourVol", "Box", 2, 2, 4 );
		yourVol.setMother( myVol );
		yourVol.setPosition( 0, -5, 0 );
		yourVol.setRotation("xyz", 0, -Math.toRadians(90), 0 );
		//yourVol.setRotation("xyz", 0, 0, -Math.toRadians(45) );
		
		Geant4Basic theirVol = new Geant4Basic("theirVol", "Box", 2, 2, 4 );
		theirVol.setMother( yourVol );
		theirVol.setPosition( 5, 0, 0 );
		theirVol.setRotation("xyz", 0, 0, -Math.toRadians(45) );
		
		System.out.println( Util.toString( topVol ) );
		
		Geant4Basic myOtherVol = Util.cloneNoChildren( myVol );
		myOtherVol.setName( "myOtherVol" );
		
		Util.moveChildrenToMother( myVol );
		
		myOtherVol.setMother( topVol );
		System.out.println( Util.toString( topVol ) );*/
		
		/*Geant4Basic originVol = new Geant4Basic("origin", "Orb", 0.02 );
		originVol.setMother( topVol );
		Geant4Basic v0Vol = Util.createArrow("v0", Util.toVector3D( v0 ), 0.2, 0.1, false, true, false );
		v0Vol.setMother( topVol );
		Geant4Basic v1Vol = Util.createArrow("v1", Util.toVector3D( v1 ), 0.2, 0.1, false, true, false );
		v1Vol.setMother( topVol );*/
		
		/*IGdmlExporter gdmlTest = VolumeExporterFactory.createGdmlFactory();
		//gdmlTest.setVerbose( true ); // not useful for large numbers of volumes
		gdmlTest.setPositionLoc("local");
		gdmlTest.setRotationLoc("local");
		gdmlTest.addTopVolume( topVol );
		gdmlTest.writeFile("test_matrix");
		System.exit(0);*/
		
		/*PrintStream originalStream = System.out; // http://stackoverflow.com/questions/8363493/hiding-system-out-print-calls-of-a-class
		PrintStream dummyStream    =  new PrintStream( new OutputStream(){
		    public void write(int b) {
		        //NO-OP
		    }} );
					
		Point3D originPos = new Point3D( 0.0, 1.0, 0.0 );
		Plane3D fidPln0 = new Plane3D( new Point3D( originPos ), new Vector3D( 0.0, 1.0, 0.0 ) );
		Plane3D fidPln1 = new Plane3D( new Point3D( fidPln0.point().x() + 1.0, fidPln0.point().y(), fidPln0.point().z() ), new Vector3D( 0.0, 0.0, 1.0 ) );
		//double[] translationShift = Util.toDoubleArray( fidPln0.point().vectorTo( fidPln1.point() ) );
		//double[] rotationShift = Util.convertVectorDiffToAxisAngle( fidPln0.normal(), fidPln1.normal() );
		
		System.out.printf(" O: %5.1f %5.1f %5.1f\n", originPos.x(), originPos.y(), originPos.z() );
		System.out.printf("F0: %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f\n", fidPln0.point().x(), fidPln0.point().y(), fidPln0.point().z(), fidPln0.normal().x(), fidPln0.normal().y(), fidPln0.normal().z() );
		System.out.printf("F1: %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f\n", fidPln1.point().x(), fidPln1.point().y(), fidPln1.point().z(), fidPln1.normal().x(), fidPln1.normal().y(), fidPln1.normal().z() );
		//System.out.printf(" S: %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f | %5.1f\n", translationShift[0], translationShift[1], translationShift[2], rotationShift[0], rotationShift[1], rotationShift[2], Math.toDegrees(rotationShift[3]) );
		
		Point3D arbPos0 = new Point3D( originPos.x(), originPos.y() + 1.0, originPos.z() + 1.0 );
		Point3D arbPos1 = new Point3D( arbPos0 );
		
		Vector3D translationVec = Util.toVector3D( translationShift );
		Vector3D centerVec = fidPln0.point().toVector3D();
		Vector3D rotationVec = new Vector3D( rotationShift[0], rotationShift[1], rotationShift[2] );
		System.out.printf("A0: %5.1f %5.1f %5.1f\n", arbPos0.x(), arbPos0.y(), arbPos0.z() );
		
		centerVec.scale( -1 );
		arbPos1.set( arbPos1, centerVec );
		System.out.printf("A1: %5.1f %5.1f %5.1f centered\n", arbPos1.x(), arbPos1.y(), arbPos1.z() );
		
		System.setOut(dummyStream); // suppress unwanted debug output from Vector3D.rotate()
		rotationVec.rotate( arbPos1, rotationShift[3] );
		System.setOut(originalStream);
		System.out.printf("A1: %5.1f %5.1f %5.1f rotated\n", arbPos1.x(), arbPos1.y(), arbPos1.z() );
		
		centerVec.scale( -1 );
		arbPos1.set( arbPos1, centerVec );
		System.out.printf("A1: %5.1f %5.1f %5.1f centered back\n", arbPos1.x(), arbPos1.y(), arbPos1.z() );
		
		arbPos1.set( arbPos1, translationVec );
		System.out.printf("A1: %5.1f %5.1f %5.1f translated\n", arbPos1.x(), arbPos1.y(), arbPos1.z() );
		
		
		Geant4Basic topVol = new Geant4Basic("top", "Box", 0 );
		
		double orbR = 0.02, arrowR = 0.01;
		
		Geant4Basic zeroVol = new Geant4Basic("origin", "Orb", orbR );
		zeroVol.setMother( topVol );
		
		Geant4Basic originVol = new Geant4Basic("origin", "Orb", orbR );
		originVol.setPosition( originPos.x()*0.1, originPos.y()*0.1, originPos.z()*0.1 );
		originVol.setMother( topVol );
		
		Geant4Basic fidPln0Vol = Util.createArrow("fidPln0", fidPln0.normal(), orbR*10, arrowR*10, true, true, false );
		fidPln0Vol.setPosition( fidPln0.point().x()*0.1, fidPln0.point().y()*0.1, fidPln0.point().z()*0.1 );
		fidPln0Vol.setMother( topVol );
		
		Geant4Basic fidPln1Vol = Util.createArrow("fidPln1", fidPln1.normal(), orbR*10, arrowR*10, true, true, false );
		fidPln1Vol.setPosition( fidPln1.point().x()*0.1, fidPln1.point().y()*0.1, fidPln1.point().z()*0.1 );
		fidPln1Vol.setMother( topVol );
		
		Geant4Basic arbPos0Vol = new Geant4Basic("arbPos0", "Orb", orbR );
		arbPos0Vol.setPosition( arbPos0.x()*0.1, arbPos0.y()*0.1, arbPos0.z()*0.1 );
		arbPos0Vol.setMother( topVol );
		
		Geant4Basic arbPos1Vol = new Geant4Basic("arbPos1", "Orb", orbR );
		arbPos1Vol.setPosition( arbPos1.x()*0.1, arbPos1.y()*0.1, arbPos1.z()*0.1 );
		arbPos1Vol.setMother( topVol );
		
		IGdmlExporter gdmlTest = VolumeExporterFactory.createGdmlExporter();
		//gdmlTest.setVerbose( true ); // not useful for large numbers of volumes
		gdmlTest.setPositionLoc("local");
		gdmlTest.setRotationLoc("local");
		gdmlTest.addTopVolume( topVol );
		gdmlTest.writeFile("test_shift");
		
		System.exit(0);*/
		
		/*Vector3D vec = new Vector3D( 0.0, 1.0, 0.5 ).asUnit();
		Matrix vecMatMul = Utils.convertRotationVectorToMatrix( vec.theta(), vec.phi() );
		Matrix vecMatCalc = Matrix.convertRotationFromEulerInZYX_ExXYZ( 0.0, vec.theta(), vec.phi() );
		double[] xyzMatMul = Matrix.convertRotationToEulerInXYZ_ExZYX( vecMatMul );
		double[] xyzMatCalc = Matrix.convertRotationToEulerInXYZ_ExZYX( vecMatCalc );
		
		Point3D pos1 = new Point3D( 0.0, 0.0, 1.0 );
		Transformation3D eulerTrans = new Transformation3D().rotateZ(xyzMatCalc[2]).rotateY(xyzMatCalc[1]).rotateX(xyzMatCalc[0]); // extrinsic ZYX
		eulerTrans.apply( pos1 );
		Point3D pos2 = new Point3D( 0.0, 0.0, 1.0 );
		Transformation3D vecTrans = new Transformation3D().rotateY(vec.theta()).rotateZ(vec.phi());
		vecTrans.apply(pos2);
		Point3D pos3 = new Point3D( 0.0, 0.0, 1.0 );
		Vector3D vecAxis = new Vector3D( -1.0, 0.0, 0.0 ).asUnit();
		vecAxis.rotate( pos3, Math.toRadians(63.4)); // axis angle
		
		vec.show();
		System.out.printf("theta % 6.1f\nphi   % 6.1f", Math.toDegrees(vec.theta()), Math.toDegrees(vec.phi()) );
		System.out.println("\nmatmul");
		vecMatMul.show();
		System.out.printf("\neuler % 6.1f % 6.1f % 6.1f", Math.toDegrees(xyzMatMul[0]), Math.toDegrees(xyzMatMul[1]), Math.toDegrees(xyzMatMul[2]) );
		System.out.println("\ncalc");
		vecMatCalc.show();
		System.out.printf("\neuler % 6.1f % 6.1f % 6.1f", Math.toDegrees(xyzMatCalc[0]), Math.toDegrees(xyzMatCalc[1]), Math.toDegrees(xyzMatCalc[2]) );
		System.out.println("\npoint 1");
		pos1.show();
		System.out.println("\npoint 2");
		pos2.show();
		System.out.println("\npoint 3");
		pos3.show();
		System.exit(0);*/
		
		//double[][] nominalData = new double[][]{ new double[]{ 1,0,-1 }, new double[]{ -1,0,-1 }, new double[]{ 0,0,2 } };
		//double[] rotationShift = new double[]{ 1,0,0,90 };
		//double[] translationShift = new double[]{ 0,1,0 };
		//double[][] measuredData = new double[][]{ new double[]{ 1,2,0 }, new double[]{ -1,2,0 }, new double[]{ 0,-1,0 } };
		//double[][] uncertainData = new double[][]{ new double[]{ 0.1,0,0 }, new double[]{ 0.1,0,0 }, new double[]{ 0,-0.1,0 } };
		//double uncertainty = 0.100; // 20 um
		
		/*double[][] nominalData = new double[][]{ new double[] { -17.350, -68.283, -286.541 }, new double[]{ 17.350, -68.283, -286.541 }, new double[]{ 3.500, -68.283, 122.333 } };
		double[][] measuredData = new double[][]{ new double[]{ -17.295, -68.028, -286.365 }, new double[]{ 17.385, -67.934, -286.308 }, new double[]{ 3.942, -68.334, 122.473 } };
		
		for( int k = 0; k < nominalData.length/3; k+=3 )
		{
			System.out.println();
			for( int j = 0; j < 3; j++ )
			{
				System.out.printf("NP%d % 8.3f % 8.3f % 8.3f", j, nominalData[k+j][0], nominalData[k+j][1], nominalData[k+j][2] );
				//System.out.printf("    UP%d % 8.3f % 8.3f % 8.3f", j, measuredData[k+j][0], measuredData[k+j][1], measuredData[k+j][2] );
				//for( int i = 0; i < 3; i++ ) measuredData[k+j][i] = measuredData[k+j][i] + uncertainData[k+j][i];
				System.out.printf("    MP%d % 8.3f % 8.3f % 8.3f", j, measuredData[k+j][0], measuredData[k+j][1], measuredData[k+j][2] );
				System.out.println();
			}
		}
		
		double[][] nominalDistances = new double[nominalData.length/3][3];
		double[][] centerData = new double[nominalData.length/3][3];
		
		for( int k = 0; k < nominalData.length/3; k+=3 ) // triangle
		{
			Point3D[] nominalPos3Ds = new Point3D[3];
			
			System.out.printf("\n%3s", " ");
			for( int j = 0; j < 3+1; j++ ) // point, with offset
			{
				if( j < 3 )	nominalPos3Ds[j] =  new Point3D(nominalData[k+j][0], nominalData[k+j][1], nominalData[k+j][2]);
				if( j > 0 ) // offset to wait for first point to be defined
				{
					if( k == 0) System.out.printf("%4s(%d %d)", " ", j-1, j%3 );
					nominalDistances[k][j-1] = nominalPos3Ds[j-1].distance( nominalPos3Ds[j%3] );
				}
			}
			
			Triangle3D centerTri = new Triangle3D( nominalPos3Ds[0], nominalPos3Ds[1], nominalPos3Ds[2] );
			centerData[k] = Util.toDoubleArray( centerTri.center().toVector3D() );
			
			System.out.printf("\nND%d % 8.3f % 8.3f % 8.3f", k, nominalDistances[k][0], nominalDistances[k][1], nominalDistances[k][2] );
			System.out.printf("    NC%d % 8.3f % 8.3f % 8.3f", k, centerData[k][0], centerData[k][1], centerData[k][2] );
			System.out.println();
		}
		
		AlignmentFactory.VERBOSE = true;
		//for(int i = 0; i < 1; i++) AlignmentFactory.adjustData( measuredData, uncertainty, nominalDistances );
		
		//System.exit(0);
		
		double[][] shiftData = AlignmentFactory.calcShifts( 1, nominalData, measuredData );
		//double[][] shiftData = AlignmentFactory.calcShifts( 1, nominalData, measuredData, uncertainty );
		
		double[][] shiftedData = AlignmentFactory.applyShift( nominalData, shiftData, centerData, 1, 1 );
		
		double[][] deltasData = AlignmentFactory.calcDeltas( 3, 3, measuredData, shiftedData );
		
		System.exit(0);*/
		
		// ConstantProvider cp = new DatabaseLoader.getSVTConstants();
		/*DatabaseConstantProvider cp = new DatabaseConstantProvider( 10, "default");
		cp.loadTable( SVTConstants.getCcdbPath() +"svt");
		cp.loadTable( SVTConstants.getCcdbPath() +"region");
		cp.loadTable( SVTConstants.getCcdbPath() +"support");
		cp.loadTable( SVTConstants.getCcdbPath() +"fiducial");
		//cp.loadTable( SVTConstants.getCcdbPath() +"material");
		cp.loadTable( SVTConstants.getCcdbPath() +"alignment");
		cp.disconnect();*/
		//System.out.println(cp.toString());
		//System.exit(0);
		
		SVTConstants.VERBOSE = true;
		DatabaseConstantProvider cp = SVTConstants.connect( true );
		
		/*SVTAlignmentFactory.setup( cp, "survey_ideals_reformat.dat", "survey_measured_reformat.dat" );
		double[][] dataFactoryIdeals = SVTAlignmentFactory.getFactoryIdealsFiducialData();
		
		SVTAlignmentFactory.calcShifts( dataFactoryIdeals, SVTAlignmentFactory.getDataSurveyMeasured(), "shifts_survey_measured_from_factory_ideals.dat" );
		//AlignmentFactory.VERBOSE = true;
		double[][] dataDeltas = SVTAlignmentFactory.calcDeltas( dataFactoryIdeals, SVTAlignmentFactory.getDataSurveyMeasured(), "deltas_survey_measured_from_factory_ideals.dat" );
		//AlignmentFactory.VERBOSE = false;
		
		//SVTAlignmentFactory.calcDeltas( dataNominal, SVTAlignmentFactory.getDataSurveyIdeals(), "deltas_survey_ideals_from_factory_nominal.dat");
		//SVTAlignmentFactory.calcDeltas( SVTAlignmentFactory.getDataSurveyIdeals(), SVTAlignmentFactory.getDataSurveyMeasured(), "deltas_survey_measured_from_survey_ideals.dat");
		
		SVTAlignmentFactory.calcTriangleSides( dataFactoryIdeals, 0, "sides_factory_ideals");
		SVTAlignmentFactory.calcTriangleSides( SVTAlignmentFactory.getDataSurveyMeasured(), 0.020, "sides_survey_measured");
		SVTAlignmentFactory.calcTriangleSides( dataDeltas, 0.020, "sides_survey_measured_from_factory_ideals");*/
		//System.exit(0);
		
		
		
		int regionSelector = 1, sectorSelector = 6;
		
		SVTVolumeFactory svtIdeal = new SVTVolumeFactory( cp, false );
		//svtIdeal.setRange( regionSelector, sectorSelector, sectorSelector );
		svtIdeal.setRange( regionSelector, 0, 0 );
		svtIdeal.makeVolumes();
		
		SVTStripFactory svtIdealStrips = new SVTStripFactory( cp, false );
		
		/*Geant4Basic module = svtNominal.createModule();
		module.setMother( svtNominal.getMotherVolume() );*/
		//Utils.shiftPosition( module, SVTGeant4Factory.MODULEWID/2, 0, SVTGeant4Factory.MODULELEN/2);
		
		//Geant4Basic region = svtNominal.createRegion( 0 );
		//region.setMother( svtNominal.getMotherVolume() );
		
		String fileNameIdealFiducials = "factory_fiducials_ideal.dat";
		Writer fileIdealFiducials = Util.openOutputDataFile( fileNameIdealFiducials );
		
		for( int region = svtIdeal.getRegionMin()-1; region < svtIdeal.getRegionMax(); region++ )
			for( int sector = svtIdeal.getSectorMin()[region]-1; sector < svtIdeal.getSectorMax()[region]; sector++ )
			{
				for( int module = svtIdeal.getModuleMin()-1; module < svtIdeal.getModuleMax(); module++ )
				{
					for( int strip = 0; strip < SVTConstants.NSTRIPS; strip+=32 )
					{
						Line3D stripLine = svtIdealStrips.getIdealStrip( region, sector, module, strip );
						//stripLine.show();
						Geant4Basic stripVol = Util.createArrow("strip"+strip+"_m"+module+"_s"+sector+"_r"+region, stripLine.toVector(), 0.5, 0.2, false, true, false );
						stripVol.setPosition( stripLine.origin().x()*0.1, stripLine.origin().y()*0.1, stripLine.origin().z()*0.1 );
						stripVol.setMother( svtIdeal.getMotherVolume() );
						//System.out.println( stripVol.gemcString() );
						//for( int c = 0; c < stripVol.getChildren().size(); c++ )
							//System.out.println( stripVol.getChildren().get(c).gemcString() );
					}
					
					Point3D[] layerCorners = svtIdealStrips.getIdealLayerCorners( region, sector, module );
					for( int i = 0; i < layerCorners.length; i++ )
					{
						Geant4Basic cornerBall = new Geant4Basic("cornerBall"+i+"_m"+module+"_s"+sector+"_r"+region, "Orb", 0.075 ); // cm
						cornerBall.setPosition( layerCorners[i].x()*0.1, layerCorners[i].y()*0.1, layerCorners[i].z()*0.1 ); // mm -> cm
						cornerBall.setMother( svtIdeal.getMotherVolume() );
					}
				}
				
				Point3D fidPos3Ds[] = SVTAlignmentFactory.getIdealFiducials( region, sector );
				
				for( int fid = 0; fid < SVTConstants.NFIDUCIALS; fid++ )
				{
					Util.writeLine( fileIdealFiducials, String.format("R%dS%02dF%d % 8.3f % 8.3f % 8.3f\n", region+1, sector+1, fid+1, fidPos3Ds[fid].x(), fidPos3Ds[fid].y(), fidPos3Ds[fid].z() ) );
					
					Geant4Basic fidBall = new Geant4Basic("fiducialBall"+fid+"_s"+sector+"_r"+region, "Orb", 0.2 ); // cm
					fidBall.setPosition( fidPos3Ds[fid].x()*0.1, fidPos3Ds[fid].y()*0.1, fidPos3Ds[fid].z()*0.1 ); // mm->cm
					fidBall.setMother( svtIdeal.getMotherVolume() );
				}
				
				Triangle3D fidTri3D = new Triangle3D( fidPos3Ds[0], fidPos3Ds[1], fidPos3Ds[2] );
				Vector3D fidVec3D = fidTri3D.normal().asUnit();
				fidVec3D.scale( 10 ); // length of arrow in mm
				Geant4Basic fidCen = Util.createArrow( "fiducialCenter_s"+sector+"_r"+region, fidVec3D, 2.0, 1.0, true, true, false );
				fidCen.setPosition( fidTri3D.center().x()*0.1, fidTri3D.center().y()*0.1, fidTri3D.center().z()*0.1 );
				fidCen.setMother( svtIdeal.getMotherVolume() );
			}
		
		Util.closeOutputDataFile( fileNameIdealFiducials, fileIdealFiducials );
		
		//System.out.println( svtIdeal.toString() );
		
		IGdmlExporter gdmlFile = VolumeExporterFactory.createGdmlFactory();
		//gdmlFile.setVerbose( true ); // not useful for large numbers of volumes
		gdmlFile.setPositionLoc("local");
		gdmlFile.setRotationLoc("local");
		gdmlFile.addTopVolume( svtIdeal.getMotherVolume() );
		gdmlFile.addMaterialPreset("mat_sensorActive", "mat_vacuum");
		gdmlFile.replaceAttribute( "structure", "volume", "name", "vol_sensorActive", "materialref", "ref", "mat_sensorActive");
		gdmlFile.replaceAttribute( "structure", "volume", "name", "vol_deadZone", "materialref", "ref", "mat_sensorActive");
		//gdmlFile.replaceAttribute( "structure", "volume", "name", "vol_rohacell", "materialref", "ref", "mat_rohacell");
		gdmlFile.writeFile("SVTFactory_ideal");
		
		//System.exit( 0 );
		
		
		
		SVTVolumeFactory svtShifted = new SVTVolumeFactory( cp, false );
		//SVTConstants.loadAlignmentShifts("shifts_survey_measured_from_factory_ideals.dat");
		//SVTConstants.loadAlignmentShifts("shifts_custom.dat");
		SVTConstants.loadAlignmentShifts("shifts_test.dat");
		//SVTConstants.loadAlignmentShifts("shifts_zero.dat");
		//SVTConstants.loadAlignmentShifts( cp );
		svtShifted.setApplyAlignmentShifts( true );
		//svtShifted.setAlignmentShiftScale( 1, 10 );
		
		//System.exit( 0 );
		
		//svtShifted.setRange( regionSelector, sectorSelector, sectorSelector );
		svtShifted.setRange( regionSelector, 0, 0 );
		
		//AlignmentFactory.VERBOSE = true;
		svtShifted.makeVolumes();
		//AlignmentFactory.VERBOSE = false;
		
		SVTStripFactory svtShiftedStrips = new SVTStripFactory( cp, false );
		svtShiftedStrips.setApplyAlignmentShifts( true );
		
		String fileNameShiftedFiducials = "factory_fiducials_Shifted.dat";
		Writer fileShiftedFiducials = Util.openOutputDataFile( fileNameShiftedFiducials );
		
		for( int region = svtShifted.getRegionMin()-1; region < svtShifted.getRegionMax(); region++ )
			for( int sector = svtShifted.getSectorMin()[region]-1; sector < svtShifted.getSectorMax()[region]; sector++ )
			{
				for( int module = svtShifted.getModuleMin()-1; module < svtShifted.getModuleMax(); module++ )
				{
					for( int strip = 0; strip < SVTConstants.NSTRIPS; strip+=32 )
					{
						Line3D stripLine = svtShiftedStrips.getShiftedStrip( region, sector, module, strip );
						//stripLine.show();
						Geant4Basic stripVol = Util.createArrow("strip"+strip+"_m"+module+"_s"+sector+"_r"+region, stripLine.toVector(), 0.5, 0.2, false, true, false );
						stripVol.setPosition( stripLine.origin().x()*0.1, stripLine.origin().y()*0.1, stripLine.origin().z()*0.1 );
						stripVol.setMother( svtShifted.getMotherVolume() );
						//System.out.println( stripVol.gemcString() );
						//for( int c = 0; c < stripVol.getChildren().size(); c++ )
							//System.out.println( stripVol.getChildren().get(c).gemcString() );
					}
					
					Point3D[] layerCorners = svtShiftedStrips.getShiftedLayerCorners( region, sector, module );
					for( int i = 0; i < layerCorners.length; i++ )
					{
						Geant4Basic cornerBall = new Geant4Basic("cornerBall"+i+"_m"+module+"_s"+sector+"_r"+region, "Orb", 0.075 ); // cm
						cornerBall.setPosition( layerCorners[i].x()*0.1, layerCorners[i].y()*0.1, layerCorners[i].z()*0.1 ); // mm -> cm
						cornerBall.setMother( svtShifted.getMotherVolume() );
					}
				}
				
				/*Point3D fidPos3Ds[] = SVTAlignmentFactory.getShiftedFiducials( region, sector );
				
				for( int fid = 0; fid < SVTConstants.NFIDUCIALS; fid++ )
				{
					Util.writeLine( fileShiftedFiducials, String.format("R%dS%02dF%d % 8.3f % 8.3f % 8.3f\n", region+1, sector+1, fid+1, fidPos3Ds[fid].x(), fidPos3Ds[fid].y(), fidPos3Ds[fid].z() ) );
					
					Geant4Basic fidBall = new Geant4Basic("fiducialBall"+fid+"_s"+sector+"_r"+region, "Orb", 0.2 ); // cm
					fidBall.setPosition( fidPos3Ds[fid].x()*0.1, fidPos3Ds[fid].y()*0.1, fidPos3Ds[fid].z()*0.1 ); // mm->cm
					fidBall.setMother( svtShifted.getMotherVolume() );
				}*/
				
				//System.out.println("r"+region+"s"+sector+"k"+SVTGeant4Factory.convertRegionSector2SvtIndex( region, sector ));
				
				//Point3D fidShiftedPos3Ds[] = SVTAlignmentFactory.getShiftedFiducials( region, sector );
				
				// calculate shifted fiducials manually to show intermediate steps
				Point3D[] fidIdealPos3Ds = SVTAlignmentFactory.getIdealFiducials( region, sector ); // lab frame
				Triangle3D fidIdealTri3D = new Triangle3D( fidIdealPos3Ds[0], fidIdealPos3Ds[1], fidIdealPos3Ds[2] );
				double [] shift = SVTConstants.getAlignmentShiftData()[SVTConstants.convertRegionSector2SvtIndex( region, sector )].clone();
				Point3D[] fidShiftedPos3Ds = new Point3D[SVTConstants.NFIDUCIALS];
				
				int n = 1;
				double d = shift[6]/n;
				for( int i = 1; i < n+1; i++ )
				{
					//System.out.println("fid "+ i );
					shift[6] = i*d;
										
					for( int fid = 0; fid < SVTConstants.NFIDUCIALS; fid++ )
					{
						fidShiftedPos3Ds[fid] = new Point3D( fidIdealPos3Ds[fid] );  // reset for next step
						SVTAlignmentFactory.applyShift( fidShiftedPos3Ds[fid], shift, fidIdealTri3D.center() );
						
						if( i == n )
						{							
							if( fid == SVTConstants.NFIDUCIALS-1 )
							{
								Triangle3D fidTri3D = new Triangle3D( fidShiftedPos3Ds[0], fidShiftedPos3Ds[1], fidShiftedPos3Ds[2] );
								Vector3D fidVec3D = fidTri3D.normal().asUnit();
								fidVec3D.scale( 10 );
								Geant4Basic fidCen = Util.createArrow( "fiducialCenter_s"+sector+"_r"+region, fidVec3D, 2.0, 1.0, true, true, false );
								fidCen.setPosition( fidTri3D.center().x()*0.1, fidTri3D.center().y()*0.1, fidTri3D.center().z()*0.1 );
								fidCen.setMother( svtShifted.getMotherVolume() );
							}
						}
						
						//Geant4Basic fidBall = new Geant4Basic("fiducialBall"+fid+"_s"+sector+"_r"+region+"_"+i, "Orb", 0.2 ); // cm
						Geant4Basic fidBall = new Geant4Basic("fiducialBall"+fid+"_s"+sector+"_r"+region, "Orb", 0.2 ); // cm
						fidBall.setPosition( fidShiftedPos3Ds[fid].x()*0.1, fidShiftedPos3Ds[fid].y()*0.1, fidShiftedPos3Ds[fid].z()*0.1 ); // mm -> cm
						fidBall.setMother( svtShifted.getMotherVolume() );
					}
				}
				
				Triangle3D fidTri3D = new Triangle3D( fidShiftedPos3Ds[0], fidShiftedPos3Ds[1], fidShiftedPos3Ds[2] );
				Vector3D fidVec3D = fidTri3D.normal().asUnit();
				fidVec3D.scale( 10 ); // length of arrow in mm
				Geant4Basic fidCen = Util.createArrow( "fiducialCenter_s"+sector+"_r"+region, fidVec3D, 2.0, 1.0, true, true, false );
				fidCen.setPosition( fidTri3D.center().x()*0.1, fidTri3D.center().y()*0.1, fidTri3D.center().z()*0.1 );
				fidCen.setMother( svtShifted.getMotherVolume() );
			}
		
		Util.closeOutputDataFile( fileNameShiftedFiducials, fileShiftedFiducials );
		
		//System.out.println( svtShifted.toString() );
		
		IGdmlExporter gdmlFile2 = VolumeExporterFactory.createGdmlFactory();
		gdmlFile2.setPositionLoc("local");
		gdmlFile2.setRotationLoc("local");	
		gdmlFile2.addTopVolume( svtShifted.getMotherVolume() );
		gdmlFile2.addMaterialPreset("mat_sensorActive", "mat_vacuum");
		gdmlFile2.replaceAttribute( "structure", "volume", "name", "vol_sensorActive", "materialref", "ref", "mat_sensorActive");
		gdmlFile2.replaceAttribute( "structure", "volume", "name", "vol_deadZone", "materialref", "ref", "mat_sensorActive");
		gdmlFile2.replaceAttribute( "structure", "volume", "name", "vol_rohacell", "materialref", "ref", "mat_rohacell");
		gdmlFile2.writeFile("SVTFactory_shifted");
		
		for( int region = svtShifted.getRegionMin()-1; region < svtShifted.getRegionMax(); region++ ) // SVTGeant4Factory.NREGIONS
			for( int sector = svtShifted.getSectorMin()[region]-1; sector < svtShifted.getSectorMax()[region]; sector++ ) // SVTGeant4Factory.NSECTORS[region]
			{
				//System.out.println("r"+region+"s"+sector+"k"+SVTGeant4Factory.convertRegionSector2SvtIndex( region, sector ));
				
				Point3D fidShiftedPos3Ds[] = SVTAlignmentFactory.getShiftedFiducials( region, sector );
				
				// calculate shifted fiducials manually to show intermediate steps
				//Point3D[] fidNominalPos3Ds = SVTAlignmentFactory.getNominalFiducials( region, sector ); // lab frame
				//Triangle3D fidNominalTri3D = new Triangle3D( fidNominalPos3Ds[0], fidNominalPos3Ds[1], fidNominalPos3Ds[2] );
				//double [] shift = SVTConstants.getAlignmentShiftData()[SVTConstants.convertRegionSector2SvtIndex( region, sector )].clone();
				//Point3D[] fidShiftedPos3Ds = new Point3D[SVTConstants.NFIDUCIALS];
				
				//int n = 1;
				//double d = shift[6]/n;
				//for( int i = 1; i < n+1; i++ )
				//{
					//System.out.println("fid "+ i );
					//shift[6] = i*d;
										
					for( int fid = 0; fid < SVTConstants.NFIDUCIALS; fid++ )
					{
						//fidShiftedPos3Ds[fid] = new Point3D( fidNominalPos3Ds[fid] );  // reset for next step
						//SVTAlignmentFactory.applyShift( fidShiftedPos3Ds[fid], shift, fidNominalTri3D.center() );	
						
						//if( i == n )
						//{
							
							//Vector3D fidDiffVec3D = fidShiftedPos3Ds[fid].vectorFrom( Util.toVector3D( SVTAlignmentFactory.getDataSurveyMeasured()[SVTGeant4Factory.convertRegionSectorFid2SurveyIndex(region, sector, fid)] ).toPoint3D() );
							//Util.writeLine( fileShiftedFiducialsDeltas, String.format("R%dS%02dF%d % 8.3f % 8.3f % 8.3f\n", region+1, sector+1, fid+1, fidDiffVec3D.x(), fidDiffVec3D.y(), fidDiffVec3D.z() ) );
							
							if( fid == SVTConstants.NFIDUCIALS-1 )
							{
								Triangle3D fidTri3D = new Triangle3D( fidShiftedPos3Ds[0], fidShiftedPos3Ds[1], fidShiftedPos3Ds[2] );
								Vector3D fidVec3D = fidTri3D.normal().asUnit();
								fidVec3D.scale( 10 );
								Geant4Basic fidCen = Util.createArrow( "fiducialCenter_s"+sector+"_r"+region, fidVec3D, 2.0, 1.0, true, true, false );
								fidCen.setPosition( fidTri3D.center().x()*0.1, fidTri3D.center().y()*0.1, fidTri3D.center().z()*0.1 );
								fidCen.setMother( svtShifted.getMotherVolume() );
							}
						//}
						
						//Geant4Basic fidBall = new Geant4Basic("fiducialBall"+fid+"_s"+sector+"_r"+region+"_"+i, "Orb", 0.2 ); // cm
						Geant4Basic fidBall = new Geant4Basic("fiducialBall"+fid+"_s"+sector+"_r"+region, "Orb", 0.2 ); // cm
						fidBall.setPosition( fidShiftedPos3Ds[fid].x()*0.1, fidShiftedPos3Ds[fid].y()*0.1, fidShiftedPos3Ds[fid].z()*0.1 ); // mm -> cm
						fidBall.setMother( svtShifted.getMotherVolume() );
					//}
					}
			}
		

		Geant4Basic svtMergeVol = new Geant4Basic("merge", "Box", 0 );
		svtShifted.appendName("_shifted");
		svtIdeal.getMotherVolume().setMother( svtMergeVol );
		svtShifted.getMotherVolume().setMother( svtMergeVol );
		
		IGdmlExporter gdmlFile3 = VolumeExporterFactory.createGdmlFactory();
		gdmlFile3.setPositionLoc("local");
		gdmlFile3.setRotationLoc("local");
		gdmlFile3.addTopVolume( svtMergeVol );
		//gdmlFile3.setVerbose( true );
		gdmlFile3.addMaterialPreset("mat_shifted", "mat_vacuum");
		gdmlFile3.replaceAttribute( "structure", "volume", "name", "_shifted", "materialref", "ref", "mat_shifted");
		gdmlFile3.replaceAttribute( "structure", "volume", "name", "fiducial", "materialref", "ref", "mat_shifted");
		gdmlFile3.replaceAttribute( "structure", "volume", "name", "sectorBall", "materialref", "ref", "mat_shifted");
		gdmlFile3.replaceAttribute( "structure", "volume", "name", "rotAxis", "materialref", "ref", "mat_shifted");
		gdmlFile3.writeFile("SVTFactory_merge");
		
		
		// verify shifted fiducials against survey measured using alignment algorithm
		//SVTAlignmentFactory.calcDeltas( SVTAlignmentFactory.getDataSurveyMeasured(), SVTAlignmentFactory.getShiftedFiducialData(), "deltas_factory_shifted_from_survey_measured.dat");
		
		System.out.println("done");
	}

}
