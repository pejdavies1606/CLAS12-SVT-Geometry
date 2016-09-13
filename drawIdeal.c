
#include <TColor>

void drawIdeal() {
    gSystem->Load("libGeom");
    gSystem->Load("libGdml");
    TGeoManager *geom = TGeoManager::Import("SVTFactory_ideal.gdml");

    TList *matList = geom->GetListOfMaterials();
    TIter matNext( matList );
    cout << "setting material transparencies\n";
    int transparencyAir = 0;
    int transparencySensor = 0;
    while( mat = (TGeoMaterial*) matNext() )
    {
        TString *matName = new TString( mat->GetName() );
        if( matName->Contains("sensorActive") )
        {
            cout << "sensor " << transparencySensor;
            mat->SetTransparency( transparencySensor );
        }
        else
        {
            cout << "default " << transparencyAir;
            mat->SetTransparency( transparencyAir );
        }
        cout << " " << matName->Data() << "\n";
    }
    
    TGeoVolume *top = geom->GetTopVolume();
    //top->SetLineColor( kWhite );
    geom->SetTopVisible();

    //Int_t myBlackIndex = TColor::GetFreeColorIndex();
    //TColor *myBlack = new TColor(myBlackIndex, 0.0, 0.0, 0.0);
    
    TObjArray *volList = geom->GetListOfVolumes();
    TIter volNext( volList );
    cout << "setting volume colours\n";
    //cout << "modules red\n";
    
    while( vol = (TGeoVolume*) volNext() )
    {
        vol->SetVisContainers(kTRUE);
        TString *volName = new TString( vol->GetName() );
        if( volName->Contains("module") )
        {
            //cout << "module";
            vol->SetLineColor( kRed );
        }
        else if( volName->Contains("sensorActive") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kBlue );
        }
        else if( volName->Contains("fiducial") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kCyan );
        }
        else if( volName->Contains("rohacell") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kWhite );
        }
        else if( volName->Contains("heatSink") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kOrange );
        }
        else if( volName->Contains("heatSinkCu") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kOrange+1 );
        }
        else if( volName->Contains("carbonFiber") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kGray+3 ); // kBlack
        }
        else if( volName->Contains("busCable") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kGray+2 ); // kBlack
        }
        else
        {
            //cout << "default";
            vol->SetLineColor( kGray );
        }
        //cout << " " << volName->Data() << "\n";
    }

    cout << "VisLevel " << geom->GetVisLevel() << "\n";
    cout << "VisOption " << geom->GetVisOption() << "\n";
    //geom->SetVisLevel( 3 );
    //geom->SetVisOption( 1 );
    //cout << "new VisLevel " << geom->GetVisLevel() << "\n";
    //cout << "new VisOption " << geom->GetVisOption() << "\n";
    top->Draw("ogl");
}
