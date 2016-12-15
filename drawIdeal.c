
#include <TColor>

void drawIdeal() {
    gSystem->Load("libGeom");
    gSystem->Load("libGdml");
    TGeoManager *geom = TGeoManager::Import("SVTFactory_ideal.gdml");

    TList *matList = geom->GetListOfMaterials();
    TIter matNext( matList );
    cout << "setting material transparencies\n";
    int transparencyHide = 100;
    int transparencyHalf = 100;
    int transparencyShow = 0;
    
    while( mat = (TGeoMaterial*) matNext() )
    {
        TString *matName = new TString( mat->GetName() );
        if( matName->Contains("hide") )
        {
            cout << "hide " << transparencyHide;
            mat->SetTransparency( transparencyHide );
        }
        else if( matName->Contains("half") )
        {
            cout << "half " << transparencyHalf;
            mat->SetTransparency( transparencyHalf );
        }
        else
        {
            cout << "show " << transparencyShow;
            mat->SetTransparency( transparencyShow );
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
            vol->SetLineColor( kBlue-1 );
        }
        else if( volName->Contains("sensorPhysical") )
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
        else if( volName->Contains("carbonFiberCu") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kGray+3 );
        }
        else if( volName->Contains("carbonFiberPk") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kGray+3 );
        }
        else if( volName->Contains("busCable") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kGray+2 );
        }
        else if( volName->Contains("pcBoard") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kOrange-4 );
        }
        else if( volName->Contains("chip") )
        {
            //cout << "sensorActive";
            vol->SetLineColor( kCyan+1 );
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
