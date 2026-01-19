#include "ICRP110PhantomConstruction.hh"
#include "ICRP110PhantomConstruction.hh"
#include "ICRP110PhantomNestedParameterisation.hh"
#include "ICRP110PhantomMaterial_Female.hh"
#include "G4SystemOfUnits.hh"
#include "ICRP110PhantomMaterial_Male.hh"
#include <array>
#include <string>
#include <iostream>
#include <filesystem>
#include <vector>
// Total count: 53 materials
std::array<std::string, 53> materialNames = {
    "air",               // [00] 
    "teeth",             // [01] 
    "bone",              // [02] 
    "humeri_upper",      // [03] 
    "humeri_lower",      // [04] 
    "arm_lower",         // [05] 
    "hand",              // [06] 
    "clavicle",          // [07] 
    "cranium",           // [08] 
    "femora_upper",      // [09] 
    "femora_lower",      // [10] 
    "leg_lower",         // [11] 
    "foot",              // [12] 
    "mandible",          // [13] 
    "pelvis",            // [14] 
    "ribs",              // [15] 
    "scapulae",          // [16] 
    "spine_cervical",    // [17] 
    "spine_thoratic",    // [18] 
    "spine_lumbar",      // [19] 
    "sacrum",            // [20] 
    "sternum",           // [21] 
    "hf_upper",          // [22] 
    "hf_lower",          // [23] 
    "med_lowerarm",      // [24] 
    "med_lowerleg",      // [25] 
    "cartilage",         // [26] 
    "skin",              // [27] 
    "blood",             // [28] 
    "muscle",            // [29] 
    "liver",             // [30] 
    "pancreas",          // [31] 
    "brain",             // [32] 
    "heart",             // [33] 
    "eye",               // [34] 
    "kidney",            // [35] 
    "stomach",           // [36] 
    "intestine_sml",     // [37] 
    "intestine_lrg",     // [38] 
    "spleen",            // [39] 
    "thyroid",           // [40] 
    "bladder",           // [41] 
    "ovaries_testes",    // [42] 
    "adrenals",          // [43] 
    "oesophagus",        // [44] 
    "misc",              // [45] 
    "uterus_prostate",   // [46] 
    "lymph",             // [47] 
    "breast_glandular",  // [48] 
    "breast_adipose",    // [49] 
    "lung",              // [50] 
    "gastro_content",    // [51] 
    "urine"              // [52] 
};
std::string PenelopeMatToMCGPU( const std::string mcgpuFilePath, 
                                const std::string materialFilePath){
    std::stringstream ss;
    ss<<5000<<" "<<120005<<std::endl;
    ss<<23002<<std::endl;
    ss<<materialFilePath<<std::endl;
    ss<<mcgpuFilePath<<std::endl;

    

    std::string fullCommand = "printf \"" + ss.str() + "\" | ./create_material_debug_rodrigo.x";
    
    std::cout<<std::endl<<std::endl;
    std::cout<<fullCommand.c_str();
    std::cout<<std::endl<<std::endl;
    int result = std::system(fullCommand.c_str());
    if (result != 0){
        std::cerr << "ERROR: Failed to create " <<mcgpuFilePath<<  std::endl;
    }
    return mcgpuFilePath;

}
std::string CreateMaterialFiles(G4Material* material,const std::string target_folder_relative)
{
  
    std::stringstream ss;
    ss<<1<<std::endl;
    ss<<material->GetName()<<std::endl;
    ss<<material->GetNumberOfElements()<<std::endl;
    ss<<2<<std::endl;
    const G4ElementVector* elements = material->GetElementVector();
    const G4double* massFractions = material->GetFractionVector();
    for(size_t i=0; i<material->GetNumberOfElements(); ++i)
    {
        const G4Element* element = (*elements)[i];
        G4double fraction = massFractions[i];
        ss<<element->GetZ()<<" "<<fraction<<std::endl;
    }
    ss<<2<<std::endl;
    ss<<material->GetDensity()/ (g/cm3)<<std::endl;
    ss<<2<<std::endl;
    std::string materialFilePath = material->GetName() + ".mat";
    ss<<materialFilePath<<std::endl;

    // std::cout<<std::endl<<std::endl;
    // std::cout<<ss.str();
    // std::cout<<std::endl<<std::endl;


    
    std::string fullCommand = "printf \"" + ss.str() + "\" | ./material.x";
    int result = std::system(fullCommand.c_str());
    if (result != 0){
        std::cerr << "ERROR: Failed to create " <<materialFilePath<<  std::endl;
    }
    ss.str("");
    ss.clear();

    std::string mcgpuFilePath = target_folder_relative + "icrp_" + material->GetName() + ".mcgpu";
    return PenelopeMatToMCGPU(mcgpuFilePath,materialFilePath);

}
void G4MaterialVectorToMCGPUFiles(const std::vector<G4Material*>& materials,const std::string target_folder_relative)
{
    //creo la lista per collegare Id e file
    std::vector<std::string> materialFileIdList;
    materialFileIdList.push_back("material/air__5-120keV.mcgpu density=0.0012 voxelId=0");
    for(size_t i = 1; i < materials.size(); ++i)
    {   
        G4Material* material = materials[i];
        std::cout << "Index: " << i << " | Material: " << material->GetName() << std::endl;

        std::string outFileName;
        if(material) {
            outFileName=CreateMaterialFiles(material,target_folder_relative);
        }
        materialFileIdList.push_back("material/"+outFileName+" "+
            "density="+std::to_string(material->GetDensity()/ (g/cm3))+" voxelId="+std::to_string(i+1));
    }

    std::string outputFileName = target_folder_relative+"MC-GPU_material_config.txt";
    std::ofstream outFile(outputFileName);

    if (outFile.is_open()) {
        for (const auto& entry : materialFileIdList) {
            outFile << entry << "\n";
        }
        outFile.close();
        std::cout << "Successfully wrote material list to " << outputFileName << std::endl;
    }
    else {
        std::cerr << "Error: Could not open file for writing!" << std::endl;
    }
}
std::vector<G4Material*> GetMaterialsMale()
{
    auto fMaterial_Male = new ICRP110PhantomMaterial_Male();
    fMaterial_Male -> DefineMaterials();

    std::vector<G4Material*> materials;
    for(size_t i = 1; i < materialNames.size(); ++i)
    {
        const std::string& name = materialNames[i];
        G4Material* material = fMaterial_Male->GetMaterial(name);
        if(material) 
            materials.push_back(material);
    }
    return materials;
}
std::vector<G4Material*> GetMaterialsFemale()
{
    auto fMaterial_Female = new ICRP110PhantomMaterial_Female();
    fMaterial_Female -> DefineMaterials();

    std::vector<G4Material*> materials;
    for(size_t i = 1; i < materialNames.size(); ++i)
    {
        const std::string& name = materialNames[i];
        G4Material* material = fMaterial_Female->GetMaterial(name);
        if(material) 
            materials.push_back(material);
    }
    return materials;
}
int main()
{
    //cartella pendb, partendo dalla cartella build di questo progetto
    //i due eseguibili dentro pendb sono material.x (per creare il .mat di penelope)
    // e create_material_debug_rodrigo.x per creare il file .mcgpu
    std::string penelopeDir = "../../../penelope/pendbase/";

    //mi sposto in pendb
    std::system("ls");
    std::filesystem::current_path(penelopeDir);
    std::system("ls");

    G4MaterialVectorToMCGPUFiles(GetMaterialsFemale(),"G4_ICRP_female/");
    G4MaterialVectorToMCGPUFiles(GetMaterialsMale(),"G4_ICRP_male/");

    return 0;
}