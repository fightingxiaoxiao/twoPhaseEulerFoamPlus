

const fvPatchList &patches = mesh.boundary();
forAll(patches, i)
{
    const std::string name = static_cast<std::string>(patches[i].name()); // 强制转换std::string
    std::smatch match;
    std::regex e("(.snow)"); // 匹配包含".snow"字符串的边界
    if (std::regex_search(name, match, e))
    {
        Info << "Recognized snow surface: \"" << patches[i].name() << "\"." << endl;
        const labelUList &faceCells = patches[i].faceCells();

        const tmp<volSymmTensorField> &Reff1 = phase1.turbulence().devRhoReff();
        const tmp<volSymmTensorField> &Reff2 = phase2.turbulence().devRhoReff();

        const vectorField &Sfp = mesh.Sf().boundaryField()[i];       // 面积矢量
        const scalarField &magSfp = mesh.magSf().boundaryField()[i]; // 面积矢量模长

        const symmTensorField &Reffp1 = Reff1().boundaryField()[i];
        const symmTensorField &Reffp2 = Reff2().boundaryField()[i];

        scalarField &alpha1p = alpha1.boundaryFieldRef()[i];
        scalarField &alpha2p = alpha2.boundaryFieldRef()[i];

        tmp<vectorField> ssp = (-Sfp / magSfp) & (alpha1p * Reffp1 + alpha2p * Reffp2);

        tmp<scalarField> uShear = sqrt(mag(ssp()) / fluid.rho()->boundaryField()[i]);

        Info << average(fluid.rho()->boundaryField()[i]) << endl;

        forAll(faceCells, i)
        {
            
            scalar deltaAlpha;
            //Info << uShear()[i] << endl;
            if (uShear()[i] > 0.27) // 侵蚀
            {
                deltaAlpha = 7e-4 * (sqr(uShear()[i]) - sqr(0.27)) * magSfp[i] / mesh.V()[faceCells[i]] * runTime.deltaTValue() / 917.;
            }
            else // 沉积
            {
                deltaAlpha = 0.;
            }

            alpha1[faceCells[i]] += deltaAlpha;
            alpha1p[i] += deltaAlpha;

            alpha2[faceCells[i]] -= deltaAlpha;
            alpha2p[i] -= deltaAlpha;

            //Info << alpha1[faceCells[i]] << endl;
        }
    }
}
