#ifndef ClusterCovMat_h
#define ClusterCovMat_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "UTIL/LCRelationNavigator.h"
#include <EVENT/MCParticle.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "DDMarlinCED.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TMatrixD.h"

using namespace lcio ;
using namespace marlin ;

class ClusterCovMat : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new ClusterCovMat;
		}
		ClusterCovMat();
		virtual ~ClusterCovMat() = default;
		ClusterCovMat(const ClusterCovMat&) = delete;
		ClusterCovMat& operator=(const ClusterCovMat&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void getNeutralCovMatFirstOrder( ReconstructedParticleImpl* NeutralPFO , std::vector<float> CovMat );
		virtual void getNeutralCovMatSecondOrder( ReconstructedParticleImpl* NeutralPFO , std::vector<float> CovMat );
		std::vector<float> getNeutralCovMatFirstOrder( TVector3 clusterPosition , float pfoEc , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();

	private:
		std::string				m_inputPfoCollection{};
		std::string				m_outputPfoCollection{};

		bool					m_AssumeNeutralPFOMassive = true;
		bool					m_isClusterEnergyKinEnergy = false;
		bool					m_updatePFO4Momentum = false;
		bool					m_useTrueJacobian = false;

};
#endif
