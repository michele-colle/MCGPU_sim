#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PenelopeRayleighModel.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4ForceCondition.hh"
#include <iostream>
#include <vector>
#include <iomanip>

int main() {
// 1. Setup the material (e.g., Aluminum)
    // 1. Setup Material
    G4Material* material = nullptr;
    try
    {
        material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    
    
    std::cout << "Material: " << material->GetName() 
              << ", Density: " << material->GetDensity()/g*cm3 << " g/cm3"
              << ", Number of Elements: " << material->GetNumberOfElements() << std::endl;

    // 2. Initialize Model
    G4PenelopeRayleighModel* penModel = new G4PenelopeRayleighModel();
    G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
    G4DataVector dummy;
    
    // Set high energy limit to ensure model is active
    penModel->SetHighEnergyLimit(150*keV);
    penModel->Initialise(gamma, dummy);

    // 3. CRITICAL: Penelope model loads atomic data only when a Cross Section is requested
    penModel->CrossSectionPerVolume(material, gamma, 5*keV);


    // 4. Trigger the built-in Geant4 dump
    // This method internally calls BuildFormFactorTable(material) if it's missing.
    std::cout << "\n--- STARTING GEANT4 PENELOPE DUMP ---" << std::endl;
    auto matName = material->GetName();
    std::replace(matName.begin(), matName.end(), '/', '_');
    std::replace(matName.begin(), matName.end(), ' ', '_');

    std::ofstream out("RayleighFF_" + matName + ".dat");
    auto oldBuf = G4cout.rdbuf();
    G4cout.rdbuf(out.rdbuf());
    penModel->DumpFormFactorTable(material);
    G4cout.rdbuf(oldBuf);
    out.close();
    std::cout << "--- DUMP COMPLETE ---" << std::endl;

    // 2. The grid setup
    int numPoints = 128;
    double xMax = 22.0; 
    
    std::vector<double> x_grid(numPoints);
    std::vector<double> f2_vals(numPoints);
    std::vector<double> p_vals(numPoints, 0.0);

    std::cout << std::scientific << std::setprecision(10);
    //std::cout << "# SECTION: PENELOPE RAYLEIGH SAMPLING TABLE (Z=" << Z << ")" << std::endl;

    for (int i = 0; i < numPoints; ++i) {
        // Linear grid in sqrt(x) to match MC-GPU/PENELOPE spacing
        double sqrtX = (std::sqrt(xMax) / (numPoints - 1)) * i;
        x_grid[i] = sqrtX * sqrtX;

        // PENELOPE physics uses q in units of 1/cm or 1/Angstrom
        // q = 4*pi * sin(theta/2)/lambda
        // In this context, q (1/A) = 4 * pi * sqrt(X)
        double q = 4.0 * pi * std::sqrt(x_grid[i]);

        // 3. Get the Form Factor from the Penelope Model
        // GetSquaredFormFactor is a public method in Penelope model for a given q
        // Note: Check your G4 version; some versions use GetSquaredFormFactor(Z, q_squared)
        double F2 = 1;
        // double F2 = penRayleigh->GetSquaredFormFactor(element, q);
        f2_vals[i] = F2;

        // 4. Integrate for Cumulative Probability P
        if (i > 0) {
            double dx = x_grid[i] - x_grid[i-1];
            p_vals[i] = p_vals[i-1] + 0.5 * (f2_vals[i] + f2_vals[i-1]) * dx;
        }
    }

    // 5. Normalize and Output RITA-style
    double pTotal = p_vals.back();
    for (int i = 0; i < numPoints; ++i) {
        double P = p_vals[i] / pTotal;
        
        // Calculate A parameter (Rational Interpolation coefficient)
        // A matches the slope of the PDF at the bin boundary
        double A = 0.0;
        if (i < numPoints - 1) {
            double dx = x_grid[i+1] - x_grid[i];
            double normalized_slope = (f2_vals[i] / pTotal) * dx;
            A = 0.5 * (normalized_slope - 1.0); 
        }

        // B is often set to 0 or a very small curvature correction in RITA
        double B = 0.0; 
        int itl = i + 1;
        int itu = std::min(numPoints, i + 2);

        // std::cout << std::setw(18) << x_grid[i] 
        //           << std::setw(18) << P 
        //           << std::setw(18) << A 
        //           << std::setw(18) << B 
        //           << std::setw(6)  << itl 
        //           << std::setw(6)  << itu << std::endl;
    }

    return 0;
}