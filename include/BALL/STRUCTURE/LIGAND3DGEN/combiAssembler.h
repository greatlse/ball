#ifndef COMBIASSEMBLER_H
#define COMBIASSEMBLER_H

#include <BALL/STRUCTURE/LIGAND3DGEN/structureAssembler.h>

class CombiAssembler
{
public:
	CombiAssembler(TemplateDatabaseManager& data, CombiLibMap* clib);
	
	~CombiAssembler();
	
	void setScaffold(RFragment& scaffold);
	void setCombiLib(CombiLibMap& clib);
	
	void writeCombinations(BALL::SDFile& handle);
	
	
private:
	
	bool _connectClashFree(BALL::Atom &at1, BALL::Atom &at2, ConnectList& connections);
	void _checkAndConnect(RAtom& acceptor, RFragment& donor);
	void _addSet(RFragment& mol);
	
	RFragment*          _work_mol;
	CombiLibMap*        _r_groups;
	std::list< RAtom* > _r_atms;
	
	MoleculeConnector  _connector;
	ConnectionResolver _cresolv;

	std::list< RFragment* >            _current_combination;
	void _combineRecur(BALL::SDFile& handle);
	
	void newSetForCurrentCombination();
};

#endif // COMBIASSEMBLER_H