#include "options.hpp"

#include <filesystem>

namespace fs = std::filesystem;

AnalysisType parseAnalysisType(const std::string &type)
{
    if (type.compare("preprocess") == 0)
        return AnalysisType::Preprocess;
    else if (type.compare("ppbayes") == 0)
        return AnalysisType::PpBayes;
    else if (type.compare("asyncppbayes") == 0)
        return AnalysisType::AsyncPpBayes;
    else if (type.compare("gauss") == 0)
        return AnalysisType::Gauss;
    else if (type.compare("asyncgauss") == 0)
        return AnalysisType::AsyncGauss;
    else
        return AnalysisType::Unknown;
}

MatrixXd parseVarianceComponentsFromString(const std::string &string)
{
    static const std::string groupSeparator = ";";
    static const std::string componentSeparator = ",";

    Gadget::Tokenizer groupTokenizer;
    groupTokenizer.getTokens(string, groupSeparator);

    // If we find no semi-colon, assume there is only one group
    if (groupTokenizer.empty())
        groupTokenizer.push_back(string);

    assert(!groupTokenizer.empty());

    Gadget::Tokenizer componentTokenizer;
    componentTokenizer.getTokens(groupTokenizer[0], componentSeparator);

    if (componentTokenizer.empty()) {
        cout << "Failed to parse variance components: " << string << endl;
        return {};
    }

    MatrixXd S(groupTokenizer.size(), componentTokenizer.size());
    const auto expectedComponentCount = static_cast<Gadget::Tokenizer::size_type>(S.cols());
    for (auto group = 0; group < S.rows(); ++group) {
        componentTokenizer.getTokens(groupTokenizer[group], componentSeparator);
        if (componentTokenizer.size() != expectedComponentCount) {
            cout << "Incorrect number of variance components! "
                 << "Got: " << componentTokenizer.size()
                 << ", expected: " << expectedComponentCount
                 << endl;
            return {};
        }

        for (unsigned int i = 0; i < expectedComponentCount; ++ i) {
            try {
                S(group, i) = stod(componentTokenizer[i]);
            }
            catch (const std::invalid_argument &) {
                cerr << "Could not parse variance component: " << componentTokenizer[i] << endl;
                return {};
            }
            catch (const std::out_of_range &) {
                cerr << "Variance component is out of range: " << componentTokenizer[i] << endl;
                return {};
            }
        }
    }

    return S;
}

MatrixXd parseVarianceComponentsFromFile(const fs::path &path)
{
    ifstream in(path);
    if (!in.is_open()) {
        cout << "Error opening variance components file: " << path << endl;
        return {};
    }

    return parseVarianceComponentsFromString({istreambuf_iterator<char>(in), istreambuf_iterator<char>()});
}

MatrixXd Options::parseVarianceComponents(const std::string &arg)
{
    const fs::path path(arg);
    if (fs::is_regular_file(path))
        return parseVarianceComponentsFromFile(path);
    else
        return parseVarianceComponentsFromString(arg);
}

void Options::inputOptions(const int argc, const char* argv[]){
    stringstream ss;
    for (unsigned i=1; i<argc; ++i) {
        if (!strcmp(argv[i], "--inp-file")) {
            optionFile = argv[++i];
            readFile(optionFile);
            return;
        } else {
            if (i==1) ss << "\nOptions:\n\n";
        }
        if (!strcmp(argv[i], "--analysis-type")) {
            const auto type = argv[++i];
            analysisType = parseAnalysisType(type);

            ss << "--analysis-type" << type << "\n";
        }
        else if (!strcmp(argv[i], "--preprocess")) {
            analysisType = AnalysisType::Preprocess;
            ss << "--preprocess " << "\n";
        }
        else if (!strcmp(argv[i], "--compress")) {
            compress = true;
            ss << "--compress " << "\n";
        }
        else if (!strcmp(argv[i], "--data-file")) {
            dataFile = argv[++i];
            inputType = getInputType(dataFile);
            ss << "--data-file " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--pheno")) {
            phenotypeFile = argv[++i];
            ss << "--pheno " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--mcmc-samples")) {
            mcmcSampleFile = argv[++i];
            ss << "--mcmc-samples " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--chain-length")) {
            chainLength = atoi(argv[++i]);
            ss << "--chain-length " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--burn-in")) {
            burnin = atoi(argv[++i]);
            ss << "--burn-in " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--seed")) {
            seed = atoi(argv[++i]);
            ss << "--seed " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--out")) {
            title = argv[++i];
            ss << "--out " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--thin")) {
            thin = atoi(argv[++i]);
            ss << "--thin " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--S")) {
            S = parseVarianceComponents(argv[++i]);
            ss << "--S " << argv[i] << "\n";
        }
        //Daniel group assignment file
        else if (!strcmp(argv[i], "--group")) {
            groupFile = argv[++i];
            ss << "--group " << argv[i] << "\n";
        }
        // Fixed effects matrix file
        else if (!strcmp(argv[i], "--fixed_effects")) {
        	fixedFile = argv[++i];
           	ss << "--fixed_effects " << argv[i] << "\n";
      	}
        // Fixed effects number
        else if (!strcmp(argv[i], "--fixedEffectNumber")) {
                fixedEffectNumber = atoi(argv[++i]);
                ss << "--fixedEffectNumber " << argv[i] << "\n";
        }
        // Failure vector file
        else if (!strcmp(argv[i], "--failure")) {
        	failureFile = argv[++i];
        	ss << "--failure " << argv[i] << "\n";
		}
        else if (!strcmp(argv[i], "--bayesW_version")) {
        		bayesW_version = argv[++i];
               	ss << "--bayesW_version " << argv[i] << "\n";
       	}
        else if (!strcmp(argv[i], "--quad_points")) {
        		quad_points = argv[++i];
                ss << "--quad_points " << argv[i] << "\n";
          	}


        else if (!strcmp(argv[i], "--thread")) {
            numThread = atoi(argv[++i]);
            ss << "--thread " << argv[i] << "\n";
        }
        else if (!strcmp(argv[i], "--sparse-data")) {
            string sparseDataType = argv[++i];
            if (sparseDataType == "eigen")
                preprocessDataType = PreprocessDataType::SparseEigen;
            else if (sparseDataType == "ragged")
                preprocessDataType = PreprocessDataType::SparseRagged;
            else
                preprocessDataType = PreprocessDataType::None;

            ss << "--sparse-data " << sparseDataType << "\n";
        }
        else if(!strcmp(argv[i], "--thread-spawned")) {
            numThreadSpawned = atoi(argv[++i]);
            ss << "--thread-spawned " << argv[i] << "\n";
        }
        else if(!strcmp(argv[i], "--decompression-concurrency")) {
            decompressionNodeConcurrency = atoi(argv[++i]);
            ss << "--decompression-concurrency " << argv[i] << "\n";
        }
        else if(!strcmp(argv[i], "--decompression-tokens")) {
            decompressionTokens = atoi(argv[++i]);
            ss << "--decompression-tokens " << argv[i] << "\n";
        }
        else if(!strcmp(argv[i], "--analysis-concurrency")) {
            analysisNodeConcurrency = atoi(argv[++i]);
            ss << "--analysis-concurrency " << argv[i] << "\n";
        }
        else if(!strcmp(argv[i], "--analysis-tokens")) {
            analysisTokens = atoi(argv[++i]);
            ss << "--analysis-tokens " << argv[i] << "\n";
        }
        else if(!strcmp(argv[i], "--preprocess-chunks")) {
            preprocessChunks = atoi(argv[++i]);
            ss << "--preprocess-chunks " << argv[i] << "\n";
        }
	else if (!strcmp(argv[i], "--iterLog")) {
	    iterLog=true;
            iterLogFile = argv[++i];
            ss << "--iterLog " << argv[i] << "\n";
        }
	else if (!strcmp(argv[i], "--colLog")) {
	     colLog=true;
	     colLogFile = argv[++i];
	     ss << "--colLog " << argv[i] << "\n";
	}
        else if(!strcmp(argv[i], "--marker-cache")) {
            useMarkerCache = true;
            ss << "--marker-cache\n";
        }
        else if(!strcmp(argv[i], "--v0E")){
	    v0E = static_cast<double>(atof(argv[++i]));
	    ss << "--v0E" << argv[i] << "\n";
	}
	  else if(!strcmp(argv[i], "--s02E")){
	    s02E = static_cast<double>(atof(argv[++i]));
	    ss << "--s02E" << argv[i] << "\n";
	}
	  else if(!strcmp(argv[i], "--v0G")){
	    v0G = static_cast<double>(atof(argv[++i]));
	    ss << "--v0G" << argv[i] << "\n";
	}
	  else if(!strcmp(argv[i], "--s02G")){
	   s02G = static_cast<double>(atof(argv[++i]));
	    ss << "--s02G" << argv[i] << "\n";
	}
        else {
            stringstream errmsg;
            errmsg << "\nError: invalid option \"" << argv[i] << "\".\n";
            throw (errmsg.str());
        }
    }

    ss << "\nInput type: ";
    switch (inputType) {
    case InputType::Unknown:
        ss << "Unknown\n";
        break;

    case InputType::BED:
        ss << "BED\n";
        break;

    case InputType::CSV:
        ss << "CSV\n";
        break;
    }

    cout << ss.str() << endl;
}

void Options::readFile(const string &file){  // input options from file
    optionFile = file;
    stringstream ss;
    ss << "\nOptions:\n\n";
    ss << boost::format("%20s %-1s %-20s\n") %"optionFile" %":" %file;
    makeTitle();

    ifstream in(file.c_str());
    if (!in) throw ("Error: can not open the file [" + file + "] to read.");

    string key, value;
    while (in >> key >> value) {
        if (key == "phenotypeFile") {
            phenotypeFile = value;
        } else if (key == "analysisType") {
            analysisType = parseAnalysisType(value);
        } else if (key == "mcmcSampleFile") {
            mcmcSampleFile = value;
        } else if (key == "chainLength") {
            chainLength = stoi(value);
        } else if (key == "burnin") {
            burnin = stoi(value);
        } else if (key == "seed") {
            seed = stoi(value);
        } else if (key == "thin") {
            thin = stoi(value);
        } else if (key == "S") {
            S = parseVarianceComponents(value);
        } else if (key == "numThread") {
            numThread = stoi(value);
        } else if (key.substr(0,2) == "//" ||
                key.substr(0,1) == "#") {
            continue;
        } else {
            throw("\nError: invalid option " + key + " " + value + "\n");
        }
        ss << boost::format("%20s %-1s %-20s\n") %key %":" %value;
    }
    in.close();

    cout << ss.str() << endl;



}

void Options::makeTitle(void){
    title = optionFile;
    size_t pos = optionFile.rfind('.');
    if (pos != string::npos) {
        title = optionFile.substr(0,pos);
    }
}
