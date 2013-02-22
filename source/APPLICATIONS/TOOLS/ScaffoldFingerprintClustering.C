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

#include "version.h"

#include <locale>

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>


using namespace BALL;
using namespace boost;
using namespace std;


void simpleDistributionAnalyze(const unordered_map<long, vector<unsigned int> >& scfp_clusters_tmp, const unsigned int n_mols)
{
	unsigned int c_1 = 0;
	unsigned int c_10 = 0;
	unsigned int c_100 = 0;
	unsigned int c_1000 = 0;
	unsigned int c_1000X = 0;
	
	set<unsigned int> sizes;
	
	pair<long, vector<unsigned int> > grp;
	BOOST_FOREACH(grp, scfp_clusters_tmp)
	{
		unsigned int size = grp.second.size();
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
	cerr << "# Molecules:\t" << n_mols << endl;
	cerr << "# ScaffoldFPs:\t" << scfp_clusters_tmp.size() << endl;
	cerr << endl;
	cerr << "# Singletons:\t" << c_1 << endl;
	cerr << "# 1 - 10:\t" << c_10 << endl;
	cerr << "# 10 - 100:\t" <<  c_100 << endl;
	cerr << "# 100 - 1000:\t" <<  c_1000 << endl;
	cerr << "# > 1000:\t" <<  c_1000X << endl;
	cerr << endl;
	
	set<unsigned int>::reverse_iterator rit = sizes.rbegin();
	
	cerr << "# TOP 1\t" <<  *rit++ << endl;
	cerr << "# TOP 2\t" <<  *rit++ << endl;
	cerr << "# TOP 3\t" <<  *rit++ << endl;
	cerr << "# TOP 4\t" <<  *rit++ << endl;
	cerr << "# TOP 5\t" <<  *rit++ << endl;
	cerr << endl;
}


long generateFingerprintHash(const String& fp)
{
	locale loc;
	const collate<char>& coll = use_facet<collate<char> >(loc);
	
	return coll.hash(fp.data(), fp.data() + fp.length());
}


int main(int argc, char* argv[])
{
	CommandlineParser parpars("ScaffoldFingerprintClustering", "cluster molecules according to their scaffold fingerprints", VERSION, String(__DATE__), "Chemoinformatics");
	
	parpars.registerParameter("t", "Target library input file", INFILE, true);
	parpars.registerParameter("fp_col", "Column number for comma separated smiles input which contains the fingerprint", INT, true, -1);
	parpars.registerParameter("id_col", "Column number for comma separated smiles input which contains the molecule identifier", INT, true, -1);
	parpars.registerParameter("tc", "Tanimoto cutoff [default: 0.7]", DOUBLE, false, 0.7);
	parpars.registerParameter("l", "Number of fingerprints to read", INT, false, "0");
	parpars.registerParameter("nt", "Number of parallel threads to use. To use all possible threads enter <max> [default: 1]", STRING, false, "1");
	parpars.setSupportedFormats("t","smi, smi.gz, csv, csv.gz, txt, txt.gz, sdf, sdf.gz");
	
	String man = "";
	
	parpars.parse(argc, argv);
	
	float sim_cutoff = parpars.get("tc").toFloat();
	String infile = parpars.get("t");
	
	unsigned int limit = parpars.get("l").toInt();
	if (limit == 0)
	{
		limit = Limits<unsigned int>::max();
	}
	
	unsigned int fp_col;
	if (parpars.get("fp_col") != "-1")
	{
		fp_col = parpars.get("fp_col").toInt() - 1;
	}
	else
	{
		fp_col = -1;
	}
	
	unsigned int id_col;
	if (parpars.get("id_col") != "-1")
	{
		id_col = parpars.get("id_col").toInt() - 1;
	}
	else
	{
		id_col = -1;
	}
	
	unsigned int n_threads = 1;
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
	
	
	// ------------------------------------------------------------------------------------------
	// Create Scaffold Fingerprint clusters
	
	Log.level(10) << "++ --------------------------------------------------------" << endl;
	Log.level(10) << "++ STEP 1: Create Scaffold Fingerprint clusters" << endl;
	
	File zipped(infile, File::MODE_IN | File::MODE_BINARY);
	
	iostreams::filtering_istream gzip_in;
	gzip_in.push(iostreams::gzip_decompressor());
	gzip_in.push(zipped);
	
	String cid;
	String scfp;
	
	long scfp_hash;
	unsigned int n_mols = 0;
	
	vector<String> cids;
	vector<unsigned short> tmp_features;
	vector<vector<unsigned short> > mol_features;
	unordered_map<long, vector<unsigned int> > scfp_clusters_tmp;
	
	
	String line;
	vector<String> split;
	getline(gzip_in, line);
	while (getline(gzip_in, line))
	{
		if (!line.hasPrefix("#"))
		{
			line.split(split);
			
			cid = split[id_col].trim();
			scfp = split[fp_col].trim();
			
			cids.push_back(cid);
			
			BinaryFingerprintMethods::parseBinaryFingerprint(scfp, tmp_features, 1);
			mol_features.push_back(tmp_features);
			
			scfp_hash = generateFingerprintHash(scfp);
			
			if (scfp_clusters_tmp.find(scfp_hash)==scfp_clusters_tmp.end())
			{
				scfp_clusters_tmp.insert(make_pair(scfp_hash, vector<unsigned int>(1, n_mols)));
			}
			else
			{
				scfp_clusters_tmp[scfp_hash].push_back(n_mols);
			}
			
			++n_mols;
			
			if (n_mols >= limit)
			{
				break;
			}
		}
	}
	
	zipped.close();
	
	simpleDistributionAnalyze(scfp_clusters_tmp, n_mols);
	
	pair<long, vector<unsigned int> > grp;
	pair<unsigned int, vector<unsigned int> > cls;
	map<unsigned int, vector<unsigned int> > scfp_clusters;
	BOOST_FOREACH(grp, scfp_clusters_tmp)
	{
		scfp_clusters.insert(make_pair(grp.second[0], grp.second));
	}
	scfp_clusters_tmp.clear();
	
	/*
	map<unsigned int, vector<long> > grp_sizes;
	BOOST_FOREACH(grp, scfp_clusters_tmp)
	{
		if (grp_sizes.find(grp.second.size()) == grp_sizes.end())
		{
			grp_sizes.insert(make_pair(grp.second.size(), vector<long>(1, grp.first)));
		}
		else
		{
			grp_sizes[grp.second.size()].push_back(grp.first);
		}
	}
	
	fstream o("scfp_max.csv", ios::out);
	o << "ID cluster" << endl;
	map<unsigned int, vector<long> >::reverse_iterator size_it = grp_sizes.rbegin();
	for (unsigned int i=0; i!=10; ++i)
	{
		cerr << size_it->second.size() << endl;
		vector<unsigned int> scfp_group = scfp_clusters_tmp[size_it->second[0]];
		
		for (unsigned int j=0; j!=scfp_group.size(); ++j)
		{
			o << cids[scfp_group[j]] << " " << scfp_group.size() << endl;
		}
		
		++size_it;
		if (size_it == grp_sizes.rend())
		{
			break;
		}
	}
	o.close();
	*/
	
	// ------------------------------------------------------------------------------------------
	// Connected components decomposition
	
	Log.level(10) << "\n++ --------------------------------------------------------" << endl;
	Log.level(10) << "++ STEP 2: Reduce singletons" << endl;
	
	Options options;
	options.setDefaultInteger(BinaryFingerprintMethods::Option::BLOCKSIZE, 820);
	options.setDefaultReal(BinaryFingerprintMethods::Option::SIM_CUTOFF, sim_cutoff);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::N_THREADS, n_threads);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::MAX_CLUSTERS, 1000);
	options.setDefaultInteger(BinaryFingerprintMethods::Option::VERBOSITY, 15);
	
	BinaryFingerprintMethods bfm(options, mol_features);
	
	vector<unsigned int> m_indices;
	vector<vector<unsigned int> > ccs;
	multimap<unsigned int, unsigned int> cc_sizes;
	vector<vector<pair<unsigned int, float> > > nn_data;
	
	BOOST_FOREACH(cls, scfp_clusters)
	{
		m_indices.push_back(cls.second[0]);
	}
	sort(m_indices.begin(), m_indices.end());
	
	map<unsigned int, vector<unsigned int> > cluster_selection;
	map<unsigned int, vector<unsigned int> >::iterator cl_iter;
	
	bool success = bfm.averageLinkageMerger(m_indices, cluster_selection, 15);
	
	unsigned int min_index;
	vector<unsigned int> scfp_tmp;
	for (cl_iter=cluster_selection.begin(); cl_iter!=cluster_selection.end(); ++cl_iter)
	{
		min_index = Limits<unsigned int>::max();
		for (unsigned int i=0; i!=cl_iter->second.size(); ++i)
		{
			scfp_tmp = scfp_clusters[m_indices[cl_iter->second[i]]];
			
			if (scfp_tmp[0] < min_index)
			{
				min_index = scfp_tmp[0];
			}
		}
		
		for (unsigned int i=0; i!=cl_iter->second.size(); ++i)
		{
			if (scfp_clusters[m_indices[cl_iter->second[i]]][0]!=min_index)
			{
				scfp_clusters[min_index].insert(scfp_clusters[min_index].end(), 
								scfp_clusters[m_indices[cl_iter->second[i]]].begin(), scfp_clusters[m_indices[cl_iter->second[i]]].end());
				scfp_clusters.erase(scfp_clusters[m_indices[cl_iter->second[i]]][0]);
			}
		}
	}
	
	unsigned int sum = 0;
	unsigned int max = 0;
	unsigned int n_singletons = 0;
	unsigned int cluster_id = 0;
	cerr << scfp_clusters.size() << endl;
	
	fstream out("scfp_clustered.csv", ios::out);
	out << "ID cluster" << endl;
	
	BOOST_FOREACH(cls, scfp_clusters)
	{
		sum += cls.second.size();
		
		if (cls.second.size() == 1)
		{
			++n_singletons;
		}
		
		if (cls.second.size() > max)
		{
			max = cls.second.size();
		}
		
		for (unsigned int i=0; i!=cls.second.size(); ++i)
		{
			out << cids[cls.second[i]] << " " << cluster_id << endl;
		}
		
		++cluster_id;
	}
	out.close();
	
	cerr << endl;
	cerr << "A " << sum << endl;
	cerr << "S " << n_singletons << endl;
	cerr << "M " << max << endl;
	
	/*
	fstream out("scfp_cc.csv", ios::out);
	out << "ID cluster" << endl;
	unsigned int cluster_id = 0;
	unsigned int n_singletons =0;
	for (cl_iter=cluster_selection.begin(); cl_iter!=cluster_selection.end(); ++cl_iter)
	{
		if (cl_iter->second.size() > 1)
		{
			for (unsigned int i=0; i!=cl_iter->second.size(); ++i)
			{
				out << cids[m_indices[cl_iter->second[i]]] << " " << cluster_id << endl;
			}
		}
		else
		{
			++n_singletons;
		}
		
		++cluster_id;
	}
	out.close();
	
	cerr << "# Singletons:\t" << n_singletons << endl;
	
	/*
	cerr << ccs.size() << endl;
	cerr << ccs[ccs.size()-1].size() << endl;
	cerr << ccs[ccs.size()-2].size() << endl;
	cerr << ccs[ccs.size()-3].size() << endl;
	cerr << ccs[ccs.size()-4].size() << endl;
	cerr << ccs[ccs.size()-5].size() << endl;
	
	fstream out("scfp_cc.csv", ios::out);
	out << "ID cluster" << endl;
	unsigned int cluster_id = 0;
	for (unsigned int i=ccs.size() - 1; i!=0; --i)
	{
		if (ccs[i].size() > 1)
		{
			for (unsigned int j=0; j!=ccs[i].size(); ++j)
			{
				out << cids[m_indices[ccs[i][j]]] << " " << cluster_id << endl;
			}
		}
		
		++cluster_id;
	}
	
	out.close();
	*/
	
	return 0;
}































