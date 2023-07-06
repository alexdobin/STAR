#include "Parameters.h"
#include "ErrorWarning.h"

void Parameters::samAttributes(){//everything related to SAM attributes
    //outSAMattributes
    outSAMattrPresent.NH=false;//TODO re-write as class with constructor?
    outSAMattrPresent.HI=false;
    outSAMattrPresent.AS=false;
    outSAMattrPresent.NM=false;
    outSAMattrPresent.MD=false;
    outSAMattrPresent.nM=false;
    outSAMattrPresent.jM=false;
    outSAMattrPresent.jI=false;
    outSAMattrPresent.RG=false;
    outSAMattrPresent.MC=false;
    outSAMattrPresent.XS=false;
    outSAMattrPresent.ch=false;

    outSAMattrPresent.vA=false;
    outSAMattrPresent.vG=false;
    outSAMattrPresent.vW=false;
    outSAMattrPresent.rB=false;
    outSAMattrPresent.ha=false;

    outSAMattrPresent.CR=false;
    outSAMattrPresent.CY=false;
    outSAMattrPresent.UR=false;
    outSAMattrPresent.UY=false;
    outSAMattrPresent.CB=false;
    outSAMattrPresent.UB=false;
    outSAMattrPresent.GX=false;
    outSAMattrPresent.GN=false;
    outSAMattrPresent.gx=false;
    outSAMattrPresent.gn=false;    
    outSAMattrPresent.sM=false;
    outSAMattrPresent.sS=false;    
    outSAMattrPresent.sQ=false;
    outSAMattrPresent.sF=false;

    outSAMattrPresent.cN=false;
    
    //for quant SAM output only NH and HI flags
    outSAMattrPresentQuant=outSAMattrPresent;
    outSAMattrPresentQuant.NH=true;
    outSAMattrPresentQuant.HI=true;
    outSAMattrOrderQuant.push_back(ATTR_NH);
    outSAMattrOrderQuant.push_back(ATTR_HI);

    vector<string> vAttr1;
    if (outSAMattributes.at(0)=="None") {
    } else if (outSAMattributes.at(0)=="All"){
        vAttr1={"NH","HI","AS","nM","NM","MD","jM","jI","MC","ch"};
    } else if (outSAMattributes.at(0)=="Standard"){
        vAttr1={"NH","HI","AS","nM"};
    } else {
        vAttr1=outSAMattributes;
    };

    for (uint ii=0;ii<vAttr1.size();ii++) {
        if        (vAttr1.at(ii)== "NH") {
            outSAMattrOrder.push_back(ATTR_NH);
            outSAMattrPresent.NH=true;
        } else if (vAttr1.at(ii)== "HI") {
            outSAMattrOrder.push_back(ATTR_HI);
            outSAMattrPresent.HI=true;
        } else if (vAttr1.at(ii)== "AS") {
            outSAMattrOrder.push_back(ATTR_AS);
            outSAMattrPresent.AS=true;
        } else if (vAttr1.at(ii)== "NM") {
            outSAMattrOrder.push_back(ATTR_NM);
            outSAMattrPresent.NM=true;
        } else if (vAttr1.at(ii)== "MD") {
            outSAMattrOrder.push_back(ATTR_MD);
            outSAMattrPresent.MD=true;
        } else if (vAttr1.at(ii)== "nM") {
            outSAMattrOrder.push_back(ATTR_nM);
            outSAMattrPresent.nM=true;
        } else if (vAttr1.at(ii)== "jM") {
            outSAMattrOrder.push_back(ATTR_jM);
            outSAMattrPresent.jM=true;
        } else if (vAttr1.at(ii)== "jI") {
            outSAMattrOrder.push_back(ATTR_jI);
            outSAMattrPresent.jI=true;
        } else if (vAttr1.at(ii)== "vA") {
            outSAMattrOrder.push_back(ATTR_vA);
            outSAMattrPresent.vA=true;
        } else if (vAttr1.at(ii)== "vG") {
            outSAMattrOrder.push_back(ATTR_vG);
            outSAMattrPresent.vG=true;
        } else if (vAttr1.at(ii)== "vW") {
            outSAMattrOrder.push_back(ATTR_vW);
            outSAMattrPresent.vW=true;
        } else if (vAttr1.at(ii)== "ha") {
            outSAMattrOrder.push_back(ATTR_ha);
            outSAMattrPresent.ha=true;            
        } else if (vAttr1.at(ii)== "RG") {
            outSAMattrOrder.push_back(ATTR_RG);
            outSAMattrOrderQuant.push_back(ATTR_RG);
            outSAMattrPresent.RG=true;
        } else if (vAttr1.at(ii)== "rB") {
            outSAMattrOrder.push_back(ATTR_rB);
            outSAMattrOrderQuant.push_back(ATTR_rB);
            outSAMattrPresent.rB=true;
        } else if (vAttr1.at(ii)== "ch") {
            outSAMattrOrder.push_back(ATTR_ch);
            outSAMattrOrderQuant.push_back(ATTR_ch);
            outSAMattrPresent.ch=true;
        } else if (vAttr1.at(ii)== "MC") {
            outSAMattrOrder.push_back(ATTR_MC);
            outSAMattrOrderQuant.push_back(ATTR_MC);
            outSAMattrPresent.MC=true;
        } else if (vAttr1.at(ii)== "CR") {
            outSAMattrOrder.push_back(ATTR_CR);
            outSAMattrOrderQuant.push_back(ATTR_CR);
            outSAMattrPresent.CR=true;
        } else if (vAttr1.at(ii)== "CY") {
            outSAMattrOrder.push_back(ATTR_CY);
            outSAMattrOrderQuant.push_back(ATTR_CY);
            outSAMattrPresent.CY=true;
        } else if (vAttr1.at(ii)== "UR") {
            outSAMattrOrder.push_back(ATTR_UR);
            outSAMattrOrderQuant.push_back(ATTR_UR);
            outSAMattrPresent.UR=true;
        } else if (vAttr1.at(ii)== "UY") {
            outSAMattrOrder.push_back(ATTR_UY);
            outSAMattrOrderQuant.push_back(ATTR_UY);
            outSAMattrPresent.UY=true;
        } else if (vAttr1.at(ii)== "CB") {
            outSAMattrOrder.push_back(ATTR_CB);
            outSAMattrOrderQuant.push_back(ATTR_CB);
            outSAMattrPresent.CB=true;
        } else if (vAttr1.at(ii)== "UB") {
            outSAMattrOrder.push_back(ATTR_UB);
            outSAMattrOrderQuant.push_back(ATTR_UB);
            outSAMattrPresent.UB=true;
        } else if (vAttr1.at(ii)== "GX") {
            outSAMattrOrder.push_back(ATTR_GX);
            outSAMattrOrderQuant.push_back(ATTR_GX);
            outSAMattrPresent.GX=true;            
        } else if (vAttr1.at(ii)== "GN") {
            outSAMattrOrder.push_back(ATTR_GN);
            outSAMattrOrderQuant.push_back(ATTR_GN);
            outSAMattrPresent.GN=true;
        } else if (vAttr1.at(ii)== "gx") {
            outSAMattrOrder.push_back(ATTR_gx);
            outSAMattrOrderQuant.push_back(ATTR_gx);
            outSAMattrPresent.gx=true;            
        } else if (vAttr1.at(ii)== "gn") {
            outSAMattrOrder.push_back(ATTR_gn);
            outSAMattrOrderQuant.push_back(ATTR_gn);
            outSAMattrPresent.gn=true;            
        } else if (vAttr1.at(ii)== "sM") {
            outSAMattrOrder.push_back(ATTR_sM);
            outSAMattrOrderQuant.push_back(ATTR_sM);
            outSAMattrPresent.sM=true;
        } else if (vAttr1.at(ii)== "sS") {
            outSAMattrOrder.push_back(ATTR_sS);
            outSAMattrOrderQuant.push_back(ATTR_sS);
            outSAMattrPresent.sS=true;  
        } else if (vAttr1.at(ii)== "sF") {
            outSAMattrOrder.push_back(ATTR_sF);
            outSAMattrOrderQuant.push_back(ATTR_sF);
            outSAMattrPresent.sF=true;  
        } else if (vAttr1.at(ii)== "sQ") {
            outSAMattrOrder.push_back(ATTR_sQ);
            outSAMattrOrderQuant.push_back(ATTR_sQ);
            outSAMattrPresent.sQ=true;
        } else if (vAttr1.at(ii)== "cN") {
            outSAMattrOrder.push_back(ATTR_cN);
            outSAMattrOrderQuant.push_back(ATTR_cN);
            outSAMattrPresent.cN=true;             
        } else if (vAttr1.at(ii)== "XS") {
            outSAMattrOrder.push_back(ATTR_XS);
            outSAMattrPresent.XS=true;
            if (outSAMstrandField.type!=1) {
                inOut->logMain << "WARNING --outSAMattributes contains XS, therefore STAR will use --outSAMstrandField intronMotif" <<endl;
                outSAMstrandField.type=1;
            };
        } else {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: unknown/unimplemented SAM atrribute (tag): "<<vAttr1.at(ii) <<"\n";
            errOut <<"SOLUTION: re-run STAR with --outSAMattributes that contains only implemented attributes\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
        };
    };

    if  (!var.yes && (outSAMattrPresent.vA | outSAMattrPresent.vG)) {
        ostringstream errOut;
        errOut <<"EXITING because of fatal PARAMETER error: --outSAMattributes contains vA and/or vG tag(s), but --varVCFfile is not set\n";
        errOut <<"SOLUTION: re-run STAR with a --varVCFfile option, or without vA/vG tags in --outSAMattributes\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (!wasp.yes && outSAMattrPresent.vW) {
        ostringstream errOut;
        errOut <<"EXITING because of fatal PARAMETER error: --outSAMattributes contains vW tag, but --waspOutputMode is not set\n";
        errOut <<"SOLUTION: re-run STAR with a --waspOutputMode option, or without vW tags in --outSAMattributes\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (outSAMattrRGline[0]!="-" && !outSAMattrPresent.RG) {
        outSAMattrOrder.push_back(ATTR_RG);
        outSAMattrOrderQuant.push_back(ATTR_RG);
        outSAMattrPresent.RG=true;
        inOut->logMain << "WARNING --outSAMattrRG defines a read group, therefore STAR will output RG attribute" <<endl;
    } else if (outSAMattrRG.size()==0 && outSAMattrPresent.RG) {
            ostringstream errOut;
            errOut <<"EXITING because of FATAL INPUT ERROR: --outSAMattributes contains RG tag, but --outSAMattrRGline is not set\n";
            errOut <<"SOLUTION: re-run STAR with a valid read group parameter --outSAMattrRGline\n";
            exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };

    if (outSAMstrandField.type==1 && !outSAMattrPresent.XS) {
        outSAMattrOrder.push_back(ATTR_XS);
        inOut->logMain << "WARNING --outSAMstrandField=intronMotif, therefore STAR will output XS attribute" <<endl;
    };

    if (wasp.yes && !outSAMattrPresent.vW) {
        outSAMattrOrder.push_back(ATTR_vW);
        outSAMattrOrderQuant.push_back(ATTR_vW);
        outSAMattrPresent.vW=true;
        inOut->logMain << "WARNING --waspOutputMode is set, therefore STAR will output vW attribute" <<endl;
    };

    //TODO: BAM-only flag should all have values > 100
    samAttrRequiresBAM(outSAMattrPresent.ch, "ch");
    samAttrRequiresBAM(outSAMattrPresent.CR, "CR");
    samAttrRequiresBAM(outSAMattrPresent.CY, "CY");
    samAttrRequiresBAM(outSAMattrPresent.UR, "UR");
    samAttrRequiresBAM(outSAMattrPresent.UY, "UY");
    samAttrRequiresBAM(outSAMattrPresent.CB, "CB");
    samAttrRequiresBAM(outSAMattrPresent.UB, "UB");
    samAttrRequiresBAM(outSAMattrPresent.sM, "sM");
    samAttrRequiresBAM(outSAMattrPresent.sS, "sS");
    samAttrRequiresBAM(outSAMattrPresent.sS, "sF");
    samAttrRequiresBAM(outSAMattrPresent.sQ, "sQ");
    samAttrRequiresBAM(outSAMattrPresent.rB, "rB");
    samAttrRequiresBAM(outSAMattrPresent.vG, "vG");
    samAttrRequiresBAM(outSAMattrPresent.vA, "vA");
    samAttrRequiresBAM(outSAMattrPresent.vW, "vW");
    samAttrRequiresBAM(outSAMattrPresent.GX, "GX");
    samAttrRequiresBAM(outSAMattrPresent.GN, "GN");
    
    if (outSAMattrPresent.GX || outSAMattrPresent.GN) {//in case GX/GN were requested without Solo
        quant.gene.yes = true;
        quant.yes = true;
        pSolo.samAttrFeature = SoloFeatureTypes::Gene;//TODO: this is a hack... need to perform read annotations independent of Solo
    };
};

void Parameters:: samAttrRequiresBAM(bool attrYes, string attrTag) {
    if (!attrYes)
        return; // nothing to do
    if (!outBAMunsorted && !outBAMcoord) {
        ostringstream errOut;
        errOut <<"EXITING because of fatal PARAMETER error: --outSAMattributes contains "<< attrTag <<" tag, which requires BAM output.\n";
        errOut <<"SOLUTION: re-run STAR with --outSAMtype BAM Unsorted (and/or) SortedByCoordinate option, or without "<< attrTag <<" tag in --outSAMattributes\n";
        exitWithError(errOut.str(), std::cerr, inOut->logMain, EXIT_CODE_PARAMETER, *this);
    };
    
    if (outSAMbool)
        inOut->logMain << "WARNING: --outSAMattributes contains "<< attrTag <<" tag. It will be output into BAM file(s), but not SAM file.\n";
};