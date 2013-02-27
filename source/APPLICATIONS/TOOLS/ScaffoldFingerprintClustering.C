// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//


#include <BALL/COMMON/limits.h>
#include <BALL/DATATYPE/string.h>
#include <BALL/FORMAT/commandlineParser.h>
#include <BALL/FORMAT/lineBasedFile.h>
#include <BALL/FORMAT/molFileFactory.h>
#include <BALL/STRUCTURE/binaryFingerprintMethods.h>
#include <BALL/SYSTEM/sysinfo.h>
#include <BALL/SYSTEM/timer.h>

#include "version.h"

#include <locale>

#include <boost/unordered_map.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>


using namespace BALL;
using namespace boost;
using namespace std;


// Number of threads to use
unsigned int n_threads;

// Limit of molecules to be read
unsigned int limit;

// Mmaximum size of a final cluster
unsigned int c_max;

// Column number which contains the fingerprint
int fp_col;

// Column number which contains the fingerprint
int scfp_col;

// Column number which contains a identifier of the compounds
int id_col;


struct Molecule
{
	String m_name;
	
	unsigned int m_index;
	
	bool is_medoid;
	
	float average_sim;
	
	vector<unsigned short> m_fp;
};

unsigned int n_mols;
unordered_map<unsigned int, vector<Molecule*> > molecules;


void simpleDistributionAnalysis()
{
	unsigned int c_1 = 0;
	unsigned int c_10 = 0;
	unsigned int c_100 = 0;
	unsigned int c_1000 = 0;
	unsigned int c_1000X = 0;
	
	set<unsigned int> sizes;
	
	for (unordered_map<unsigned int, vector<Molecule*> >::iterator it=molecules.begin(); it!=molecules.end(); ++it)
	{
		unsigned int size = it->second.size();
		sizes.insert(size);
		
		if (size==1)
		{
			++c_1;
		}
		else
		{
			if (size<=10)
			{
				++c_10;
			}
			else
			{
				if (size<=100)
				{
					++c_100;
				}
				else
				{
					if (size<=1000)
					{
						++c_1000;
					}
					else
					{
						++c_1000X;
					}
				}
			}
		}
	}
	
	cerr << endl;
	cerr << "++ Molecules:\t" << n_mols << endl;
	cerr << "++ ScaffoldFPs:\t" << molecules.size() << endl;
	cerr << "++" << endl;
	cerr << "++ Singletons:\t" << c_1 << endl;
	cerr << "++ 1 - 10:\t" << c_10 << endl;
	cerr << "++ 10 - 100:\t" <<  c_100 << endl;
	cerr << "++ 100 - 1000:\t" <<  c_1000 << endl;
	cerr << "++ > 1000:\t" <<  c_1000X << endl;
	cerr << "++" << endl;
	
	unsigned int count = 0;
	cerr << "++ Max sizes: ";
	for (set<unsigned int>::reverse_iterator rit=sizes.rbegin(); rit!=sizes.rend(); ++rit)
	{
		cerr << " " << *rit;
		
		if (count == 20)
		{
			break;
		}
		
		++count;
	}
	cerr << endl;
}


void readFullFingerprints(const String& infile)
{
	vector<Molecule*> tmp_molecules;
	tmp_molecules.resize(n_mols);
	
	for (unordered_map<unsigned int, vector<Molecule*> >::iterator it=molecules.begin(); it!=molecules.end(); ++it)
	{
		for (unsigned int i=0; i!=it->second.size(); ++i)
		{
			it->second[i]->m_fp.clear();
			tmp_molecules[it->second[i]->m_index] = it->second[i];
		}
	}
	
	File zipped(infile, File::MODE_IN | File::MODE_BINARY);
	
	iostreams::filtering_istream gzip_in;
	gzip_in.push(iostreams::gzip_decompressor());
	gzip_in.push(zipped);
	
	unsigned int n_mols = 0;
	vector<unsigned short> tmp_features;
	
	String line;
	vector<String> split;
	getline(gzip_in, line);
	while (getline(gzip_in, line))
	{
		if (!line.hasPrefix("#"))
		{
			line.split(split);
			BinaryFingerprintMethods::parseBinaryFingerprint(split[fp_col].trim(), tmp_features, 1);
			
			tmp_molecules[n_mols]->m_fp = tmp_features;
			
			++n_mols;
			
			if (n_mols >= limit)
			{
				break;
			}
		}
		
	}
	
	zipped.close();
}


long generateFingerprintHash(const String& fp)
{
	locale loc;
	const collate<char>& coll = use_facet<collate<char> >(loc);
	
	return coll.hash(fp.data(), fp.data() + fp.length());
}


void scaffoldFingerprintMerging(const String& infile)
{
	File zipped(infile, File::MODE_IN | File::MODE_BINARY);
	
	iostreams::filtering_istream gzip_in;
	gzip_in.push(iostreams::gzip_decompressor());
	gzip_in.push(zipped);
	
	n_mols = 0;
	
	long hash;
	Molecule* mol;
	unsigned int index;
	vector<unsigned short> tmp_features;
	unordered_map<long, unsigned int> hash_index_mapping;
	
	String line, scfp;
	vector<String> split;
	getline(gzip_in, line);
	while (getline(gzip_in, line))
	{
		if (!line.hasPrefix("#"))
		{
			line.split(split);
			
			mol = new Molecule;
			mol->m_name = split[id_col].trim();
			mol->m_index = n_mols;
			mol->is_medoid = false;
			mol->average_sim = 1.0;
			
			scfp = split[scfp_col].trim();
			hash = generateFingerprintHash(scfp);
			
			if (hash_index_mapping.find(hash) == hash_index_mapping.end())
			{
				BinaryFingerprintMethods::parseBinaryFingerprint(scfp, tmp_features, 1);
				mol->m_fp = tmp_features;
				
				index = molecules.size();
				hash_index_mapping[hash] = index;
				molecules[index] = vector<Molecule*>(1, mol);
			}
			else
			{
				molecules[hash_index_mapping[hash]].push_back(mol);
			}
			
			++n_mols;
			
			if (n_mols % 1000 == 0)
			{
				cerr << "\r++ Molecules read:  " << n_mols << "                 ";
			}
			if (n_mols >= limit)
			{
				break;
			}
		}
	}
	cerr << "\r++ Molecules read:  " << n_mols << endl;
	
	zipped.close();
	
	return;
}


void scaffoldNearestNeighbourMerging()
{
	vector<unsigned int> m_indices;
	vector<vector<unsigned short>* > tmp_features;
	unordered_map<unsigned int, unsigned int> index_mapping;
	
	unsigned int count = 0;
	for (unordered_map<unsigned int, vector<Molecule*> >::iterator it=molecules.begin(); it!=molecules.end(); ++it)
	{
		m_indices.push_back(count);
		tmp_features.push_back(&it->second[0]->m_fp);
		
		index_mapping.insert(make_pair(count, it->first));
		
		++count;
	}
	
	Options options;
	options.setDefaultInteger(BinaryFingerprintMethods::Option::BLOCKSIZE, 820);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::N_THREADS, n_threads);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::MAX_CLUSTERS, 100);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::VERBOSITY, 15);
	
	BinaryFingerprintMethods bfm(options, tmp_features);
	
	map<unsigned int, vector<unsigned int> > cluster_selection;
	
	bfm.averageLinkageMerger(m_indices, cluster_selection, 5);
	
	unsigned int index, insert_index;
	for (map<unsigned int, vector<unsigned int> >::iterator it=cluster_selection.begin(); it!=cluster_selection.end(); ++it)
	{
		insert_index = index_mapping[it->second[0]];
		
		for (unsigned int i=1; i!=it->second.size(); ++i)
		{
			index = index_mapping[it->second[i]];
			
			molecules[insert_index].insert(molecules[insert_index].end(), molecules[index].begin(), molecules[index].end());
			molecules.erase(index);
		}
	}
	
	count = 0;
	unordered_map<unsigned int, vector<Molecule*> > id_remapped;
	for (unordered_map<unsigned int, vector<Molecule*> >::iterator it=molecules.begin(); it!=molecules.end(); ++it)
	{
		id_remapped[count++] = it->second;
	}
	
	molecules = id_remapped;
	id_remapped.clear();
	
	return;
}


void averageLinkageClustering(const vector<unsigned int>& to_cluster)
{
	Options options;
	options.setDefaultInteger(BinaryFingerprintMethods::Option::BLOCKSIZE, 820);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::N_THREADS, n_threads);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::MAX_CLUSTERS, 1000);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::VERBOSITY, 15);
	
	BinaryFingerprintMethods bfm(options);
	
	bool first_cl;
	unsigned int index;
	vector<Molecule*> tmp, c_tmp;
	vector<unsigned int> m_indices;
	vector<vector<unsigned short>* > tmp_features;
	
	vector<pair<unsigned int, float> > nn_data;
	map<unsigned int, vector<unsigned int> > cluster_selection;
	
	for (unsigned int i=0; i!=to_cluster.size(); ++i)
	{
		m_indices.clear();
		tmp_features.clear();
		tmp = molecules[to_cluster[i]];
		molecules.erase(to_cluster[i]);
		
		for (unsigned int j=0; j!=tmp.size(); ++j)
		{
			m_indices.push_back(j);
			tmp_features.push_back(&tmp[j]->m_fp);
		}
		
		bfm.setLibraryFeatures(tmp_features);
		bfm.averageLinkageClustering(m_indices, nn_data, cluster_selection);
		
		first_cl = true;
		for (map<unsigned int, vector<unsigned int> >::iterator cit=cluster_selection.begin(); cit!=cluster_selection.end(); ++cit)
		{
			c_tmp.clear();
			for (unsigned int j=0; j!=cit->second.size(); ++j)
			{
				c_tmp.push_back(tmp[cit->second[j]]);
			}
			
			if (first_cl)
			{
				first_cl = false;
				index = to_cluster[i];
			}
			else
			{
				index = molecules.size();
			}
			
			molecules[index] = c_tmp;
		}
	}
	
	return;
}


void calculateClusterMedoids()
{
	Options options;
	options.setDefaultInteger(BinaryFingerprintMethods::Option::BLOCKSIZE, 820);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::N_THREADS, n_threads);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::MAX_CLUSTERS, 1000);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::VERBOSITY, 0);
	
	BinaryFingerprintMethods bfm(options);
	
	unsigned int medoid_index;
	
	vector<float> avg_sims;
	vector<unsigned int> m_indices;
	vector<vector<unsigned short>* > tmp_features;
	
	for (unordered_map<unsigned int, vector<Molecule*> >::iterator it=molecules.begin(); it!=molecules.end(); ++it)
	{
		if (it->second.size() == 1)
		{
			it->second[0]->is_medoid = true;
		}
		else
		{
			m_indices.clear();
			tmp_features.clear();
			
			for (unsigned int i=0; i!=it->second.size(); ++i)
			{
				m_indices.push_back(i);
				tmp_features.push_back(&it->second[i]->m_fp);
			}
			
			bfm.setLibraryFeatures(tmp_features);
			bfm.calculateSelectionMedoid(m_indices, medoid_index, avg_sims);
			
			it->second[medoid_index]->is_medoid = true;
			
			for (unsigned int i=0; i!=it->second.size(); ++i)
			{
				it->second[i]->average_sim = avg_sims[i];
			}
		}
	}
}


int main(int argc, char* argv[])
{
	CommandlineParser parpars("ScaffoldFingerprintClustering", "cluster molecules according to their scaffold fingerprints", VERSION, String(__DATE__), "Chemoinformatics");
	
	parpars.registerParameter("t", "Target library input file", INFILE, true);
	parpars.registerParameter("o", "Result output file", OUTFILE, true);
	parpars.registerParameter("scfp_col", "Column number which contains the scaffold fingerprint", INT, true, -1);
	parpars.registerParameter("fp_col", "Column number which contains full the fingerprint", INT, true, -1);
	parpars.registerParameter("id_col", "Column number for comma separated smiles input which contains the molecule identifier", INT, true, -1);
	parpars.registerParameter("c_max", "Maximum size of a final cluster. [default 1000]", INT, false, 1000);
	parpars.registerParameter("l", "Number of fingerprints to read", INT, false, "0");
	parpars.registerParameter("nt", "Number of parallel threads to use. To use all possible threads enter <max> [default: 1]", STRING, false, "1");
	parpars.setSupportedFormats("t","smi, smi.gz, csv, csv.gz, txt, txt.gz, sdf, sdf.gz");
	
	parpars.setSupportedFormats("t","smi, smi.gz, csv, csv.gz, txt, txt.gz, sdf, sdf.gz");
	parpars.setSupportedFormats("o","smi, smi.gz, csv, csv.gz, txt, txt.gz, sdf, sdf.gz");
	
	String man = "";
	
	parpars.parse(argc, argv);
	
	String infile = parpars.get("t");
	
	limit = parpars.get("l").toInt();
	if (limit == 0)
	{
		limit = Limits<unsigned int>::max();
	}
	
	c_max = 1000;
	if (parpars.get("c_max") != "1000")
	{
		c_max = parpars.get("c_max").toInt();
	}
	
	if (parpars.get("scfp_col") != "-1")
	{
		scfp_col = parpars.get("scfp_col").toInt() - 1;
	}
	else
	{
		scfp_col = -1;
	}
	
	if (parpars.get("fp_col") != "-1")
	{
		fp_col = parpars.get("fp_col").toInt() - 1;
	}
	else
	{
		fp_col = -1;
	}
	
	if (parpars.get("id_col") != "-1")
	{
		id_col = parpars.get("id_col").toInt() - 1;
	}
	else
	{
		id_col = -1;
	}
	
	n_threads = 1;
	if (parpars.get("nt") != "1")
	{
		if (parpars.get("nt") == "max")
		{
			n_threads = SysInfo::getNumberOfProcessors();
		}
		else
		{
			if (parpars.get("nt").toInt() > SysInfo::getNumberOfProcessors())
			{
				n_threads = SysInfo::getNumberOfProcessors();
				Log.info() << "++ INFO: Specified number of threads exceeds available threads. Setting number to available threads." << endl;
			}
			else
			{
				n_threads = parpars.get("nt").toInt();
			}
		}
	}
	
	if (File::isAccessible(parpars.get("o")))
	{
		Log.error() << "-- FAILED: Specified output file already exists: " << parpars.get("o") << endl;
		Log.error() << endl;
		
		return 1;
	}
	
	
	
	Timer* t = new Timer();
	t->start();
	
	
	Log.setLevel(12);
	
	
	// ------------------------------------------------------------------------------------------
	// Create Scaffold Fingerprint clusters
	
	Log.level(10) << "++ --------------------------------------------------------" << endl;
	Log.level(10) << "++ STEP 1: Create Scaffold Fingerprint clusters" << endl;
	
	// Get scaffold fingerprint clusters
	scaffoldFingerprintMerging(infile);
	
	// Get some information about scaffold clusters
	simpleDistributionAnalysis();
	
	
	
	// ------------------------------------------------------------------------------------------
	// Reduce singletons
	
	Log.level(10) << "\n++ --------------------------------------------------------" << endl;
	Log.level(10) << "++ STEP 2: Reduce singletons" << endl;
	
	scaffoldNearestNeighbourMerging();
	
	// Get some information about scaffold clusters
	simpleDistributionAnalysis();
	
	
	
	// ------------------------------------------------------------------------------------------
	// Cluster large scaffoldFP pre-clusters
	
	Log.level(10) << "\n++ --------------------------------------------------------" << endl;
	Log.level(10) << "++ STEP 3: Cluster large scaffoldFP pre-clusters " << endl;
	
	vector<unsigned int> to_cluster;
	for (unordered_map<unsigned int, vector<Molecule*> >::iterator it=molecules.begin(); it!=molecules.end(); ++it)
	{
		if (it->second.size() > c_max)
		{
			to_cluster.push_back(it->first);
		}
	}
	
	readFullFingerprints(infile);
	
	if (!to_cluster.empty())
	{
		averageLinkageClustering(to_cluster);
		to_cluster.clear();
		
		// Get some information about scaffold clusters
		simpleDistributionAnalysis();
	}
	else
	{
		Log.level(10) << "++ Nothing to do ... " << endl;
	}
	
	
	
	// ------------------------------------------------------------------------------------------
	// Cluster Medoids
	
	Log.level(10) << "\n++ --------------------------------------------------------" << endl;
	Log.level(10) << "++ STEP 4: Calculate cluster medoids " << endl;
	
	calculateClusterMedoids();
	
	
	
	// ------------------------------------------------------------------------------------------
	// Remap ids and write clusters
	
	Log.level(10) << "\n++ --------------------------------------------------------" << endl;
	Log.level(10) << "++ STEP 5: Write final clusters " << endl;
	
	File out(parpars.get("o"), File::MODE_OUT);
	out << "# MolID:         external molecule identifier." << endl;
	out << "# ClusterID:     cluster identifier is a integer value [1-n] where n is the total number of clusters." << endl;
	out << "# ClusterMedoid: 1 = molecule is medoid of its cluster. 0 = not medoid of cluster." << endl;
	out << "#                Multiple medoids are possible due to duplicate fingerprints. All duplicates of a medoid are also marked as medoids." << endl;
	out << "# AverageSim:    Average similarity of fingerprint to all others in cluster." << endl;
	out << "MolID\tClusterID\tClusterMedoid\tAverageSim" << endl;
	
	for (unordered_map<unsigned int, vector<Molecule*> >::iterator it=molecules.begin(); it!=molecules.end(); ++it)
	{
		for (unsigned int i=0; i!=it->second.size(); ++i)
		{
			out << it->second[i]->m_name << "\t" << it->first << "\t" << it->second[i]->is_medoid << "\t" << it->second[i]->average_sim << endl;
			
			delete it->second[i];
		}
		it->second.clear();
	}
	molecules.clear();
	
	out.close();
	
	
	t->stop();
	LongSize seconds = t->getClockTime();
	delete t;
	
	// ------------------------------------------------------------------------------------------
	// Timing
	
	Log.level(10) << "\n++ --------------------------------------------------------" << endl;
	
	if (seconds > 60)
	{
		if (seconds > 3600)
		{
			Log.level(10) << "++ Elapsed time:  " << seconds / 3600 << " hours" << endl;
		}
		else
		{
			Log.level(10) << "++ Elapsed time:  " << seconds / 60 << " minutes" << endl;
		}
	}
	else
	{
		Log.level(10) << "++ Elapsed time:  " << seconds << " seconds" << endl;
	}
	
	
	Log.level(10) << "\n++" << endl;
	Log.level(10) << "++ DONE" << endl;
	Log.level(10) << "++" << endl;
	
	
	return 0;
}































