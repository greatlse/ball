// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
#ifndef CLASHRESOLVER_H
#define CLASHRESOLVER_H

#include <BALL/STRUCTURE/LIGAND3DGEN/ligand3DBase.h>

#include <BALL/MATHS/angle.h>
#include <BALL/DATATYPE/hashSet.h>
#include <BALL/STRUCTURE/geometricTransformations.h>


/**
 * This whole clash resolving module should be reworked to implement a fast
 * and efficient clash resolving method. Several good solutions to this problem
 * are already known, but the implementation might be at least a bachelor
 * thesis.
 */
namespace BALL 
{

/// C l a s s   C l a s h D e t e c t o r
/// ############################################################################
/**
 * @brief The ClashDetector class - Simply provides methods to detect clashes
 * within and between molecules. It simply counts all clashes, where a clash
 * is defined as two atoms having a distance that is less than their added
 * vdw-radii minus a tolerance. The standard tolerance used here is 
 * 1.2 Angsroem.
 *
 * #TODO: perhaps its smarter to define a relative tolerance here
 */
class ClashDetector
{
public:
	/**
	 * @brief ClashDetector - Initializes the ClashDetector. Here it is possible
	 * to specify an alternative tolerance value
	 *
	 * @param tolerance
	 */
	ClashDetector( float tolerance = 1.2 );
	~ClashDetector();
	
	/**
	 * @brief detectPairList - detects the number of clashes (return value)
	 * from a AtomPairList. Every pair is checked for potentially causing a clash
	 *
	 * #TODO: possibly reuse BALL::AtomBijection here!
	 *
	 * @return 
	 */
	int detectPairList(std::list< std::pair< BALL::Atom*, BALL::Atom*> >& p_list);
	
	/**
	 * @brief createBetweenPairList - creates a pairlist (p_list) for the two
	 * AtomContainer ac1 and ac2. Every atom from ac1 is checked agains every atom
	 * from ac2. The two molecules may be different. Thus we can check for clashes
	 * that occur between two molecules, but ignore clashes that exist within each
	 * molecule.
	 *
	 * @param ac1
	 * @param ac2
	 * @param p_list
	 */
	void createBetweenPairList(BALL::AtomContainer& ac1, BALL::AtomContainer& ac2,
														 std::list< std::pair< BALL::Atom*, BALL::Atom*> >& p_list);
	
	/**
	 * @brief createInnerPairList - forms a list from a single atomcontainer
	 * input. The list checks each atom of the molecule agains all remaining
	 * atoms. Thus we can check for clashes a molecule causes within itself.
	 *
	 * @param ac
	 * @param p_list
	 */
	void createInnerPairList(BALL::AtomContainer& ac, 
													 std::list< std::pair< BALL::Atom*, BALL::Atom*> >& p_list);
	
	/**
	 * @brief detectInMolecule - checks each atom of 'ac' against all 
	 * remaining atoms in 'ac' if a clash exists. Each clash is counted once and
	 * the total count is returned.
	 * 
	 * @param ac
	 * @return 
	 */
	int detectInMolecule(BALL::AtomContainer& ac);
	
	/**
	 * @brief detectBetweenMolecules - Detect all clashes between 'ac1' and 'ac2'
	 * but ignores clashes occuring within each fragment
	 * 
	 * @param ac1
	 * @param ac2
	 * @return 
	 */
	int detectBetweenMolecules(BALL::AtomContainer& ac1, BALL::AtomContainer& ac2);
	
private:
	/*
	 * @brief atom3Away - check if two atoms are speparated by at least 3 bonds
	 * 
	 * @param at1
	 * @param at2
	 * @return 
	 */
	bool atom3Away(BALL::Atom& at1, BALL::Atom& at2);
	
	/*
	 * @brief doClash - check if actual dist+tolerance is grater than the 
	 * theoretical vdw-dist
	 * 
	 * @param atm1
	 * @param atm2
	 * @return 
	 */
	bool doClash(BALL::Atom& atm1, BALL::Atom& atm2);
	
	float _tolerance;
};




/// C l a s s   R o t a t o r
/// ############################################################################
/**
 * @brief The Rotator class - The general idea is taken from Jan Fuhrmann and
 * Marcel Schumanns implementation of 'rotate' found in 
 * BALL_DOCKING_GENETICDOCK_ROTATE_BOND. The class abstracts all steps that are
 * necessary to rotate along a bond given by two atoms.
 */
class Rotator
{
public:
	Rotator();
	~Rotator();
	
	/**
	 * @brief rotate - rotate along a previously specified axis about 'angle'
	 * degree/radiant. If no axis was set before, nothing will happen.
	 * 
	 * @param angle
	 */
	void rotate(const BALL::Angle &angle);
	
	/**
	 * @brief setAxis - set the rotation axis to the positions given by 'atm1' and
	 * 'atm2' and set atoms that follow after 'atm2' to be rotated. If parent is
	 * specified, the atoms that belong to 'parent' are to be rotated.
	 * 
	 * @param atm1
	 * @param atm2
	 * @param parent
	 * @param angle
	 */
	void setAxis(BALL::Atom& atm1, BALL::Atom& atm2, BALL::Composite* parent = 0);
	
private:
	
	/*
	 * @brief init - set axis and translation from the coordinates of the input
	 * atoms
	 * 
	 * @param atm1
	 * @param atm2
	 */
	void init(BALL::Atom& atm1, BALL::Atom& atm2);
	
	/*
	 * @brief setAtomsToRotate - select the atoms that are to be rotated. First
	 * call should have 'start' and 'probe' to be the atom in the direction
	 * where we want to add atoms to '_to_rotate', whereas 'block' is the atom
	 * that marks the side which is not added.
	 * @param start
	 * @param probe
	 * @param block
	 */
	void setAtomsToRotate(BALL::Atom &start, BALL::Atom &probe, BALL::Atom &block);
	
	/*
	 * @brief findRotateDirection - check whether all atoms in 'direction' are
	 * cildren of 'parent'. If so return true.
	 * 
	 * @param block
	 * @param direction
	 * @param parent
	 * @return 
	 */
	bool findRotateDirection(BALL::HashSet<BALL::Atom *> &tested, 
													 BALL::Atom& direction, BALL::Composite* parent);
	
	BALL::HashSet<BALL::Atom *> _to_rotate; // atoms to be rotated
	
	// these define the rotation axis:
	BALL::Matrix4x4 _translation;
	BALL::Vector3 _axis;
	BALL::Vector3 _basis;
	
};




/// C l a s s   C l a s h R e s o l v e r
/// ############################################################################
/**
 * @brief The ClashResolver class - resolves clashes within a molecule
 */
class ClashResolver: public ClashDetector
{
public:
	/**
	 * @brief ClashResolver - create a clashResolver and set the tolerance
	 * for the underlying ClashDetector.
	 *
	 * @param tolerance
	 * @param max_rotors
	 */
	ClashResolver( float tolerance = 1.2);
	~ClashResolver();
	
	/**
	 * @brief setMolecule - init the ClashResolver with a molecule pointer
	 * and a list of possible rotors to solve the clashes.
	 *
	 * @param molecule
	 * @param connections
	 */
	void setMolecule(BALL::AtomContainer& molecule, BALL::ConnectList& rotors);
	
	/**
	 * @brief resolve - resolves clashes within the set molecule by rotating
	 * it along the list of rotable bonds (given via setMolecule())
	 *
	 * Here we have two times the-number-of-rotors tries to set all rotors
	 * to a value creating fewest clashes. So we apply a primitive greedy
	 * method.
	 *
	 * @return
	 */
	int resolve();
	
	/**
	 * Solve optimally under the constraint of discretized angles. ('steps'
	 * is the number by which we devide 360° to get discrete angles)
	 * 
	 * Rotates all bridging bonds
	 */ 
	int resolveAllCombinations( const int& steps = 3 );
	
	
protected:
	BALL::AtomContainer* _molecule;

	BALL::ConnectList* _all_rotors;
	
private:
	int allCombinationsRecur(const BALL::ConnectList::iterator& p, 
													const BALL::ConnectList::iterator &p_end, 
													BALL::Angle& angle, const int& steps );
};





/// C l a s s   C o n n e c t i o n R e s o l v e r
/// ############################################################################
/**
 * @brief The ConnectionResolver class - resolves clashes that occur between two
 * molecules AND clashes that were newly introduced withn a molecule
 */
class ConnectionResolver: public ClashDetector
{
public:
	
	ConnectionResolver( float tolerance = 1.2, int max_rotors = 2);
	~ConnectionResolver();
	
	/**
	 * Preconditions: 
	 * 1.) The roots of 'atm1' and 'atm2' are two different BALL::AtomContainers 
	 * 2.) A bond between 'atm1' and 'atm2' exists.
	 * 3.) 'connections' contains all (rotable) fragment-bridging bonds
	 * 3.) The AtomContainer of 'atm1' is the larger one.
	 * 4.) Optimally '_root_large' should not be transformed.
	 * @brief setMolecule, sets the Resolver up to analyze and resolve clashes
	 *        between and within the fragments behind atm1 and atm2.
	 * @param atm1 part of the larger fragment
	 * @param atm2 part of the smaller fragment
	 */
	void setMolecule(BALL::Atom& atm1, BALL::Atom& atm2, 
									 BALL::ConnectList& connection_rotors);
	
	/**
	 * @brief resolve -  does not to handle old clashes that existed already. 
	 * Always returns the correct complete clash-count, after optimization.
	 * NOTE: each resolveXYZ() call gives only clash counts for the fragment parts
	 * that should be solved but not the total clash count.
	 * 
	 * @param optimal
	 * @param conserve_large
	 * @return 
	 */
	std::pair<int, bool> resolve();
	
	/**
	 * @brief detect - Detects ONLY clashes that occur between the two fragments 
	 * that are to be connected. NOT clashes within any fragment (we assume that 
	 * given clashes can't be resolved by repeated attempts of ConnectionResolver)
	 * 
	 * @return number of clashes found BETWEEN both fragments
	 */
	int detect();
	
	static int num_resolve_calls;
	
private:
	/*
	 * Try to resolve clashes between large and small fragment by rotating along
	 * the connecting bond.
	 */
	int resolveConnection();
	
	/*
	 * Rotate bonds within fragment 'frag' to remove inter and intra clashes.
	 * 'frag' needs to be either '_large_root' or '_small_root'.
	 */
	int resolveFragment(BALL::AtomContainer& frag , BALL::ConnectList &clist, const int &given_clashes = 0);
	
	int resolveLarge(const int &given_clashes = 0);
	
	/*
	 * Solve optimally under the constraint of discretized angles. ('steps'
	 * gives the number of angles that are to be tested for each bond)
	 * 
	 * Rotates all bridging bonds in 'small' and 'large'
	 */ 
	int resolveAllCombinations( const int& steps );
	
	int allCombinationsRecur(const BALL::ConnectList::iterator& p, 
													 const BALL::ConnectList::iterator &p_end, 
													 BALL::Angle& angle, 
													 const int& steps , 
													 int &best_cnt);
	
	int detectAll();
	
	BALL::Atom* atm_large;
	BALL::Atom* atm_small;
	
	BALL::AtomContainer* _large_root;
	BALL::AtomContainer* _small_root;

	BALL::ConnectList* _small_rotors;
	BALL::ConnectList* _large_rotors;
	BALL::ConnectList* _all_rotors;
	
	BALL::TemplateCoord* _save_large;
	BALL::TemplateCoord* _save_small;
	
	std::list< std::pair<BALL::Atom*, BALL::Atom*> > _clash_small;
	std::list< std::pair<BALL::Atom*, BALL::Atom*> > _clash_large;
	std::list< std::pair<BALL::Atom*, BALL::Atom*> > _clash_inter;
	
	int _max_rotations;     // maximum number of bonds to rotate
};

}// End Namespace "BALL"
#endif // CLASHRESOLVER_H