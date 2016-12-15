
void test(){
    gSystem->Load("libGeom");
    gSystem->Load("libGdml");
    TGeoManager *geom = TGeoManager::Import("test_matrix.gdml");
    TList *matList = geom->GetListOfMaterials();
    TIter matNext( matList );
    while( mat = (TGeoMaterial*) matNext() )
    {
        mat->SetTransparency( 0 );
    }
    TObjArray *volList = geom->GetListOfVolumes();
    TIter volNext( volList );
    while( vol = (TGeoVolume*) volNext() )
    {
        vol->SetVisContainers(kTRUE);
    }
    TGeoVolume *top = geom->GetTopVolume();
    geom->SetTopVisible();
    top->Draw("ogl");
}

