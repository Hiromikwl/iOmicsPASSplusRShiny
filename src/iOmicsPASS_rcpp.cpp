/*
 * iOmicsPASS_rcpp.cpp
 *
 * Rcpp entry point for iOmicsPASS.
 * Replaces the standalone main.cpp binary call.
 *
 * Author: Hiromi WL Koh / Claude Sonnet 4.6
 */

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include "globals.hpp"
#include "geneClass.hpp"
#include "subjectClass.hpp"

#include <algorithm>
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <fstream>
#include <sstream>

using namespace std;
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Helper: read a tab-delimited file back as a Rcpp::DataFrame
// ---------------------------------------------------------------------------
static Rcpp::DataFrame read_tsv_as_df(const string& path) {
    ifstream f(path.c_str());
    if(!f.is_open()) return Rcpp::DataFrame::create();

    vector<string> header;
    string line;
    if(!getline(f, line)) return Rcpp::DataFrame::create();
    {
        istringstream ss(line);
        string tok;
        while(getline(ss, tok, '\t')) header.push_back(tok);
    }

    int ncols = (int)header.size();
    vector<vector<string>> cols(ncols);
    while(getline(f, line)){
        if(line.empty()) continue;
        istringstream ss(line);
        string tok;
        int c = 0;
        while(getline(ss, tok, '\t') && c < ncols){
            cols[c++].push_back(tok);
        }
        while(c < ncols) cols[c++].push_back("");
    }

    Rcpp::List lst(ncols);
    for(int j = 0; j < ncols; j++){
        // Try numeric conversion
        bool allnum = true;
        vector<double> num_col(cols[j].size());
        for(int i = 0; i < (int)cols[j].size(); i++){
            const string& s = cols[j][i];
            if(s == "" || s == "NA") { num_col[i] = NA_REAL; continue; }
            try { num_col[i] = stod(s); }
            catch(...) { allnum = false; break; }
        }
        if(allnum){
            Rcpp::NumericVector v(num_col.begin(), num_col.end());
            lst[j] = v;
        } else {
            Rcpp::CharacterVector v(cols[j].begin(), cols[j].end());
            lst[j] = v;
        }
    }
    lst.attr("names") = Rcpp::wrap(header);
    Rcpp::DataFrame df(lst);
    return df;
}

// ---------------------------------------------------------------------------
// Helper: read a single-column text file as a character vector
// ---------------------------------------------------------------------------
static Rcpp::CharacterVector read_lines_as_charvec(const string& path) {
    ifstream f(path.c_str());
    Rcpp::CharacterVector ret;
    if(!f.is_open()) return ret;
    string line;
    while(getline(f, line)){
        if(!line.empty()) ret.push_back(line);
    }
    return ret;
}


// ---------------------------------------------------------------------------
// iOmicsPASS_core — main Rcpp export
// ---------------------------------------------------------------------------

//' Run iOmicsPASS core algorithm from R
//'
//' @param dataX Feature x sample DataFrame for data type X (primary omics).
//'   Row 1 = feature names (col 1), remaining columns = per-sample values.
//' @param dataY Optional Feature x sample DataFrame for data type Y.
//' @param dataZ Optional Feature x sample DataFrame for data type Z (normaliser).
//' @param phenotype DataFrame with columns: SampleID, Group.
//' @param withinNet DataFrame with columns: nodeA, nodeB, direction (1/-1).
//' @param btwNet Optional between-network DataFrame: nodeA (X), nodeB (Y), direction.
//' @param pathway Optional pathway DataFrame: feature, pathway, annotation.
//' @param prior Optional prior probability DataFrame: SampleID, class columns.
//' @param log_X,log_Y,log_Z Apply log2 transform to each data type.
//' @param ztrans_X,ztrans_Y,ztrans_Z Apply z-transform to each data type.
//' @param analyzeY,analyzeZ Whether data Y / data Z were supplied.
//' @param knnimpute Use KNN imputation for missing values.
//' @param knnK K for KNN imputation.
//' @param blocksize Block size for candidate matrix in KNN.
//' @param crossV Run K-fold cross-validation.
//' @param numFold Number of CV folds.
//' @param minthres User-supplied minimum threshold (overrides CV if crossV=FALSE).
//' @param enrichment Run sub-network enrichment.
//' @param bgProp Minimum proportion of pathway genes required in background.
//' @param minbgSize Minimum number of background edges for a pathway.
//' @param minsigSize Minimum number of signature edges for enrichment output.
//' @param minobs Minimum number of non-missing observations per group per feature.
//' @param minprop Minimum proportion of non-missing observations per group per feature.
//' @param usePrior Use prior probabilities in discriminant scoring.
//' @param normalizeZby Which data (\"X\" or \"Y\") is used to normalise data Z.
//' @return Named list with analysis results.
//' @export
// [[Rcpp::export]]
Rcpp::List iOmicsPASS_core(
        Rcpp::DataFrame dataX,
        Rcpp::Nullable<Rcpp::DataFrame> dataY     = R_NilValue,
        Rcpp::Nullable<Rcpp::DataFrame> dataZ     = R_NilValue,
        Rcpp::DataFrame phenotype                 = Rcpp::DataFrame(),
        Rcpp::DataFrame withinNet                 = Rcpp::DataFrame(),
        Rcpp::Nullable<Rcpp::DataFrame> btwNet    = R_NilValue,
        Rcpp::Nullable<Rcpp::DataFrame> pathway   = R_NilValue,
        Rcpp::Nullable<Rcpp::DataFrame> prior     = R_NilValue,
        bool log_X        = false,
        bool log_Y        = false,
        bool log_Z        = false,
        bool ztrans_X     = true,
        bool ztrans_Y     = true,
        bool ztrans_Z     = true,
        bool analyzeY     = false,
        bool analyzeZ     = false,
        bool knnimpute    = false,
        int  knnK         = 10,
        int  blocksize    = 1500,
        bool crossV       = true,
        int  numFold      = 10,
        double minthres   = 1e6,
        bool enrichment   = true,
        double bgProp     = 0.5,
        int  minbgSize    = 1,
        int  minsigSize   = 3,
        int  minobs       = 1,
        double minprop    = 0.5,
        bool usePrior     = false,
        std::string normalizeZby = "Y"
) {

    // -----------------------------------------------------------------------
    // 1. Temporary output directory
    // -----------------------------------------------------------------------
    string tmpDir = Rcpp::as<string>(Rcpp::Function("tempdir")()) + "/iOmicsPASS_run/";
    Rcpp::Function dirCreate("dir.create");
    dirCreate(tmpDir,
              Rcpp::Named("showWarnings") = false,
              Rcpp::Named("recursive")    = true);

    Rcpp::Rcout << "\n\nInitiating iOmicsPASS (Rcpp)...\n" << std::endl;

    // -----------------------------------------------------------------------
    // 2. Build Input_t
    // -----------------------------------------------------------------------
    Input_t UserInput;
    UserInput.log_X        = log_X;
    UserInput.log_Y        = log_Y;
    UserInput.log_Z        = log_Z;
    UserInput.ztrans_X     = ztrans_X;
    UserInput.ztrans_Y     = ztrans_Y;
    UserInput.ztrans_Z     = ztrans_Z;
    UserInput.analyzeY     = analyzeY;
    UserInput.analyzeZ     = analyzeZ;
    UserInput.knnimpute    = knnimpute;
    UserInput.knnK         = knnK;
    UserInput.blocksize    = blocksize;
    UserInput.crossV       = crossV;
    UserInput.numFold      = numFold;
    UserInput.minthres     = minthres;
    UserInput.enrichment   = enrichment;
    UserInput.bgProp       = bgProp;
    UserInput.minbgSize    = minbgSize;
    UserInput.minsigSize   = minsigSize;
    UserInput.minobs       = minobs;
    UserInput.minprop      = minprop;
    UserInput.prior_input  = usePrior;
    UserInput.normalize_Zby = normalizeZby;
    UserInput.directory    = tmpDir;

    // -----------------------------------------------------------------------
    // 3. Phenotype → subtypeMap
    // -----------------------------------------------------------------------
    map<string, set<string> > subtypeMap;
    Rcpp::Rcout << "\nReading phenotype information..." << std::endl;
    set<string> subjects = df_to_subtypeMap(phenotype, subtypeMap);
    vector<string> grplab = getKey(&subtypeMap);

    // -----------------------------------------------------------------------
    // 4. Data → subjectMap
    // -----------------------------------------------------------------------
    map<string, subjectClass> subjectMap;
    Rcpp::Rcout << "\nReading data X..." << std::endl;
    vector<string> X_g = df_to_subjectMap(dataX, subtypeMap, subjectMap, UserInput, subjects, 3);
    Rcpp::Rcout << "There are " << X_g.size() << " features in data X." << std::endl;

    vector<string> Y_g, Z_g;
    if(analyzeY && !dataY.isNull()){
        Rcpp::Rcout << "\nReading data Y..." << std::endl;
        Y_g = df_to_subjectMap(Rcpp::as<Rcpp::DataFrame>(dataY),
                               subtypeMap, subjectMap, UserInput, subjects, 2);
        Rcpp::Rcout << "There are " << Y_g.size() << " features in data Y." << std::endl;
    }
    if(analyzeZ && !dataZ.isNull()){
        Rcpp::Rcout << "\nReading data Z..." << std::endl;
        Z_g = df_to_subjectMap(Rcpp::as<Rcpp::DataFrame>(dataZ),
                               subtypeMap, subjectMap, UserInput, subjects, 1);
        Rcpp::Rcout << "There are " << Z_g.size() << " features in data Z." << std::endl;
    }

    // -----------------------------------------------------------------------
    // 5. Network files → interaction maps
    // -----------------------------------------------------------------------
    map<string, set<string> > WithinMap, BTWMap;
    map<string, int> InteractMap;
    Rcpp::Rcout << "\nReading within-network..." << std::endl;
    set<string> WITHIN_edge = df_to_withinNet(withinNet, WithinMap, InteractMap,
                                               X_g, subjects, subjectMap);

    set<string> BTW_edge;
    vector<string> target_g;
    if(analyzeY && !btwNet.isNull()){
        if(!analyzeZ) target_g = Y_g;
        else {
            if(normalizeZby=="Y") target_g = IntersectVec2(Y_g, Z_g, 2);
            else                  X_g      = IntersectVec2(X_g, Z_g, 3);
        }
        Rcpp::Rcout << "\nReading between-network..." << std::endl;
        BTW_edge = df_to_btwNet(Rcpp::as<Rcpp::DataFrame>(btwNet), BTWMap, InteractMap,
                                target_g, X_g, subjects, subjectMap);
    } else if(!analyzeY && analyzeZ){
        X_g = IntersectVec2(X_g, Z_g, 3);
    }

    set<string> Alledges = UnionSet(WITHIN_edge, BTW_edge);

    // Write BGlist for reference
    printBGlist(tmpDir, Alledges);

    // -----------------------------------------------------------------------
    // 6. Pathways (optional)
    // -----------------------------------------------------------------------
    map<string, set<string> > PATHmap;
    map<string, string> PATHwayAnnot;
    set<string> pathwaySet;
    if(enrichment && !pathway.isNull()){
        Rcpp::Rcout << "\nReading pathway file..." << std::endl;
        pathwaySet = df_to_pathway(Rcpp::as<Rcpp::DataFrame>(pathway),
                                   PATHmap, PATHwayAnnot);
    }

    // -----------------------------------------------------------------------
    // 7. Prior probabilities (optional)
    // -----------------------------------------------------------------------
    if(usePrior && !prior.isNull()){
        Rcpp::Rcout << "\nInserting prior probabilities..." << std::endl;
        df_to_prior(Rcpp::as<Rcpp::DataFrame>(prior), subjectMap);
    }

    // -----------------------------------------------------------------------
    // 8. Determine common subjects
    // -----------------------------------------------------------------------
    vector<string> sub_X = getSubjects(3, X_g, subjectMap);
    vector<string> sub_Y, sub_Z;
    if(analyzeY) sub_Y = getSubjects(2, Y_g, subjectMap);
    if(analyzeZ) sub_Z = getSubjects(1, Z_g, subjectMap);

    vector<string> commonSub;
    if(!analyzeY) commonSub = sub_X;
    if(analyzeY)  commonSub = IntersectVec(sub_Y, sub_X);
    if(analyzeZ)  commonSub = IntersectVec(commonSub, sub_Z);

    Rcpp::Rcout << "\nNumber of common subjects: " << commonSub.size() << std::endl;

    // -----------------------------------------------------------------------
    // 9. KNN imputation (optional)
    // -----------------------------------------------------------------------
    if(knnimpute){
        bool checkP = checkNAs(subjectMap, 3);
        bool checkR = analyzeY ? checkNAs(subjectMap, 2) : false;
        bool checkD = analyzeZ ? checkNAs(subjectMap, 1) : false;

        if(checkD){
            Rcpp::Rcout << "Imputing data Z...\n";
            map<set<string>, vector<vector<double> > > CandidateSET_d;
            vector<vector<double> > MAT_D = creatematrix(Z_g, sub_Z, subjectMap, 1);
            determineCandidate(blocksize, Z_g, MAT_D, CandidateSET_d);
            KNNimpute(UserInput.cid_Z, sub_Z, Z_g, MAT_D, CandidateSET_d, subjectMap, knnK, 1);
        }
        if(checkR){
            Rcpp::Rcout << "Imputing data Y...\n";
            map<set<string>, vector<vector<double> > > CandidateSET_r;
            vector<vector<double> > MAT_R = creatematrix(Y_g, sub_Y, subjectMap, 2);
            determineCandidate(blocksize, Y_g, MAT_R, CandidateSET_r);
            KNNimpute(UserInput.cid_Y, sub_Y, Y_g, MAT_R, CandidateSET_r, subjectMap, knnK, 2);
        }
        if(checkP){
            Rcpp::Rcout << "Imputing data X...\n";
            map<set<string>, vector<vector<double> > > CandidateSET_p;
            vector<vector<double> > MAT_P = creatematrix(X_g, sub_X, subjectMap, 3);
            determineCandidate(blocksize, X_g, MAT_P, CandidateSET_p);
            KNNimpute(UserInput.cid_X, sub_X, X_g, MAT_P, CandidateSET_p, subjectMap, knnK, 3);
        }
    }

    // -----------------------------------------------------------------------
    // 10. Standardize features (z-transform)
    // -----------------------------------------------------------------------
    if(analyzeZ && ztrans_Z) StandardizeFeatures(tmpDir, commonSub, Z_g, 1, subjectMap);
    if(analyzeY && ztrans_Y) StandardizeFeatures(tmpDir, commonSub, Y_g, 2, subjectMap);
    if(ztrans_X)             StandardizeFeatures(tmpDir, commonSub, X_g, 3, subjectMap);

    // -----------------------------------------------------------------------
    // 11. Cleanup subtypeMap to only common subjects; build neighbor map
    // -----------------------------------------------------------------------
    CleanupMap(commonSub, subtypeMap);

    map<string, set<string> > Targets2BTWmap;
    if(analyzeY) ReverseMap(BTWMap, Targets2BTWmap);

    map<string, vector<set<string> > > FirstDNeighborMap;
    set<string> features = createNeighborMap(tmpDir, &X_g, &target_g,
                                              &BTWMap, &Targets2BTWmap,
                                              &WithinMap, FirstDNeighborMap);

    // -----------------------------------------------------------------------
    // 12. Fill PAM map (shrunken centroids)
    // -----------------------------------------------------------------------
    Rcpp::Rcout << "\nCarrying out shrunken centroid module..." << std::endl;
    map<string, PAM_t> PAMmap;
    bool useZ = analyzeZ;
    vector<double> thresMax = fillPAMmap2(useZ, normalizeZby, &Alledges, &commonSub,
                                           &InteractMap, &subtypeMap, &subjectMap,
                                           PAMmap, &FirstDNeighborMap, tmpDir);
    Rcpp::Rcout << "Maximum Dij (average): " << calMean(thresMax) << std::endl;

    // -----------------------------------------------------------------------
    // 13. Cross-validation
    // -----------------------------------------------------------------------
    double minThres = UserInput.minthres;
    if(crossV){
        Rcpp::Rcout << "\nCarrying out cross-validation...(This may take a while)\n" << std::endl;
        double cvThres = CVKfold(tmpDir, useZ, usePrior, normalizeZby, numFold,
                                  &Alledges, &commonSub, 30, &InteractMap,
                                  &FirstDNeighborMap, &subtypeMap, &subjectMap, &PAMmap);
        if(minThres == NA_VAL) minThres = cvThres;
    } else {
        Rcpp::Rcout << "\nSkipping cross-validation. Using threshold: " << minThres << std::endl;
    }

    // Scale per-class thresholds
    double Avgthres = calMean(thresMax);
    vector<double> MODthres = thresMax;
    for(int i = 0; i < (int)thresMax.size(); i++)
        MODthres[i] = (minThres / Avgthres) * thresMax[i];

    // -----------------------------------------------------------------------
    // 14. OutputGeneSurv — writes output files and returns edge survival sets
    // -----------------------------------------------------------------------
    vector<vector<set<string> > > geneSurv =
        outputGeneSurv(tmpDir, useZ, usePrior, normalizeZby, MODthres,
                       &Alledges, &commonSub, &subtypeMap, &PAMmap,
                       &subjectMap, &InteractMap);

    // -----------------------------------------------------------------------
    // 15. Network enrichment (optional)
    // -----------------------------------------------------------------------
    if(enrichment && !pathway.isNull()){
        Rcpp::Rcout << "\nCarrying out sub-network enrichment..." << std::endl;
        map<string, set<string> > PATHmapEdge;
        vector<string> PATHs = createEdgeOrientedPathwayMap(&features, &Alledges,
                                                              &PATHmap, &PATHmapEdge,
                                                              bgProp, minbgSize);
        set<string> EdgeSurv_grp;
        for(int j = 0; j < 2; j++){
            for(int i = 0; i < (int)grplab.size(); i++){
                EdgeSurv_grp = (geneSurv[j])[i];
                CreateOutput(tmpDir, grplab[i], PATHs, j, &EdgeSurv_grp, &Alledges,
                             &PATHmapEdge, &PATHmap, &PATHwayAnnot, minsigSize);
            }
        }
    }

    Rcpp::Rcout << "\n\nThe iOmicsPASS analysis is complete!\n" << std::endl;

    // -----------------------------------------------------------------------
    // 16. Collect results into an R List
    // -----------------------------------------------------------------------
    Rcpp::List result;

    // Background edge list
    result["BGlist"] = read_lines_as_charvec(tmpDir + "BGlist.txt");

    // Minimum threshold selected
    result["minThres"] = minThres;

    // Helper lambda to safely read a TSV or return R_NilValue
    auto safe_read = [&](const string& path) -> Rcpp::RObject {
        std::ifstream chk(path.c_str());
        if(chk.good()) return Rcpp::wrap(read_tsv_as_df(path));
        return R_NilValue;
    };

    // CV errors table
    if(crossV) result["CVerrors"] = safe_read(tmpDir + "CVerrors.txt");
    else        result["CVerrors"] = R_NilValue;

    // Selected predictive edges
    result["EdgesSelected"]      = safe_read(tmpDir + "EdgesSelected_minThres.txt");
    result["SampleProbabilities"] = safe_read(tmpDir + "SampleClass_Probabilities.txt");
    result["AttributesTable"]    = safe_read(tmpDir + "AttributesTable.txt");

    // Z-transformed data
    result["ZtransX"] = (ztrans_X)             ? safe_read(tmpDir + "Ztransform_dataX.txt") : (Rcpp::RObject)R_NilValue;
    result["ZtransY"] = (analyzeY && ztrans_Y) ? safe_read(tmpDir + "Ztransform_dataY.txt") : (Rcpp::RObject)R_NilValue;
    result["ZtransZ"] = (analyzeZ && ztrans_Z) ? safe_read(tmpDir + "Ztransform_dataZ.txt") : (Rcpp::RObject)R_NilValue;

    // Enrichment results (one per group, up/down)
    Rcpp::List enrichUp, enrichDown;
    if(enrichment && !pathway.isNull()){
        for(int i = 0; i < (int)grplab.size(); i++){
            string fUp   = tmpDir + grplab[i] + "_Enrichment_up.txt";
            string fDown = tmpDir + grplab[i] + "_Enrichment_down.txt";
            if(std::ifstream(fUp))   enrichUp[grplab[i]]   = read_tsv_as_df(fUp);
            else                     enrichUp[grplab[i]]   = R_NilValue;
            if(std::ifstream(fDown)) enrichDown[grplab[i]] = read_tsv_as_df(fDown);
            else                     enrichDown[grplab[i]] = R_NilValue;
        }
    }
    result["Enrichment_up"]   = enrichUp;
    result["Enrichment_down"] = enrichDown;

    // Store tmp directory so the user can inspect intermediate files
    result["tmpDir"] = tmpDir;

    return result;
}
