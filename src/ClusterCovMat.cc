#include "ClusterCovMat.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include <UTIL/LCRelationNavigator.h>
#include "EVENT/MCParticle.h"
#include "EVENT/Cluster.h"
#include "EVENT/ReconstructedParticle.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"
#include "marlin/VerbosityLevels.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;

ClusterCovMat aClusterCovMat;

ClusterCovMat::ClusterCovMat() :

Processor("ClusterCovMat")
{
	_description = "Set the convariance matrix in (P,E) for neutral hadrons and photons";

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"inputPfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"outputPfoCollection",
					"Name of output pfo collection",
					m_outputPfoCollection,
					std::string("updatedNeutralPFOs")
				);

	registerProcessorParameter(	"AssumeNeutralPFOMassive",
					"true: Neutral PFOs are taken massive, false: Neutral PFOs are taken massless",
					m_AssumeNeutralPFOMassive,
					bool(true)
				);
	registerProcessorParameter(	"isClusterEnergyKinEnergy",
					"true: the cluster energy is interpreted as kinetic energy of PFO, false: the cluster energy is interpreted as momentum magnitude of PFO",
					m_isClusterEnergyKinEnergy,
					bool(false)
				);

	registerProcessorParameter(	"updatePFO4Momentum",
					"true: Update 4-momentum of PFOs, false: set 4-momentum for PFOs same as input PFO",
					m_updatePFO4Momentum,
					bool(false)
				);

	registerProcessorParameter(	"useTrueJacobian",
					"true: Use (mathematically) true Jacobian for the option E_cluster = |p|, false: for the option E_cluster = |p|, Use the same jacobian as the option E_cluster = E_kinetic",
					m_useTrueJacobian,
					bool(false)
				);

}

void ClusterCovMat::init()
{

	streamlog_out(MESSAGE) << "   init called  " << std::endl;
	printParameters();

}

void ClusterCovMat::Clear()
{

}

void ClusterCovMat::processRunHeader()
{

}

void ClusterCovMat::processEvent( EVENT::LCEvent *pLCEvent )
{

	LCCollection *inputPfoCollection{};
	IMPL::LCCollectionVec* outputPfoCollection(NULL);
	outputPfoCollection = new IMPL::LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE );
	outputPfoCollection->setSubset( true );
	int n_PFO = -1;
	this->Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << pLCEvent->getEventNumber() << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		n_PFO = inputPfoCollection->getNumberOfElements();
		if ( n_PFO == -1 ) streamlog_out(DEBUG5) << "	Input PFO collection (" << m_inputPfoCollection << ") has no element (PFO) " << std::endl;
		streamlog_out(DEBUG5) << "	Total Number of PFOs: " << n_PFO << std::endl;
		for (int i_pfo = 0; i_pfo < n_PFO ; ++i_pfo)
		{
			ReconstructedParticleImpl* outputPFO = dynamic_cast<ReconstructedParticleImpl*>( inputPfoCollection->getElementAt( i_pfo ) );
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "	-------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG3) << "	Processing PFO at index " << i_pfo << std::endl;
			streamlog_out(DEBUG3) << *outputPFO << std::endl;
			float pfoMass = outputPFO->getMass();
			TVector3 clusterPosition( 0.0 , 0.0 , 0.0 );
			TLorentzVector pfoFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			double outputPFOMomentum[3]{0., 0., 0.};
			std::vector<float> outputCovMatrix( 10 , 0.0 );
			if ( ( outputPFO->getTracks() ).size() == 0 )
			{
				if ( !m_AssumeNeutralPFOMassive ) pfoMass = 0.0;
				float clusterX		= ( outputPFO->getClusters()[0] )->getPosition()[0];
				float clusterY		= ( outputPFO->getClusters()[0] )->getPosition()[1];
				float clusterZ		= ( outputPFO->getClusters()[0] )->getPosition()[2];
				clusterPosition	= TVector3( clusterX , clusterY , clusterZ );
				float clusterDistance	= sqrt( pow( clusterX , 2 ) + pow( clusterY , 2 ) + pow( clusterZ , 2 ) );
				float pfoMomentumMag	= 0;
				float pfoEnergy	= outputPFO->getEnergy();
				float pfoE;
				pfoMomentumMag = ( m_isClusterEnergyKinEnergy ? sqrt( pow( pfoEnergy , 2 ) + 2 * pfoMass * pfoEnergy ) : pfoEnergy );
				float pfoPx;
				float pfoPy;
				float pfoPz;
				if ( m_updatePFO4Momentum )
				{
					pfoPx	= pfoMomentumMag * clusterX / clusterDistance;
					pfoPy	= pfoMomentumMag * clusterY / clusterDistance;
					pfoPz	= pfoMomentumMag * clusterZ / clusterDistance;
					pfoE = ( m_isClusterEnergyKinEnergy ? pfoEnergy + pfoMass : sqrt( pow( pfoMomentumMag , 2 ) + pow( pfoMass , 2 ) ) );
				}
				else
				{
					pfoPx	= outputPFO->getMomentum()[ 0 ];
					pfoPy	= outputPFO->getMomentum()[ 1 ];
					pfoPz	= outputPFO->getMomentum()[ 2 ];
					pfoE	= outputPFO->getEnergy();
				}
				std::vector<float> clusterPositionError = ( outputPFO->getClusters()[0] )->getPositionError();
				float clusterEnergyError = ( outputPFO->getClusters()[0] )->getEnergyError();
				TVector3 pfoMomentum( pfoPx , pfoPy , pfoPz );
				pfoFourMomentum	= TLorentzVector( pfoMomentum , pfoE );
				getNeutralCovMatFirstOrder( outputPFO , outputCovMatrix );
//				outputCovMatrix	= getNeutralCovMatFirstOrder( clusterPosition , pfoEnergy , pfoMass , clusterPositionError , clusterEnergyError );

				outputPFOMomentum[ 0 ] = pfoFourMomentum.Px();
				outputPFOMomentum[ 1 ] = pfoFourMomentum.Py();
				outputPFOMomentum[ 2 ] = pfoFourMomentum.Pz();
				pfoE = pfoFourMomentum.E();


//				outputPFO->setType(outputPFO->getType());
//				outputPFO->setMomentum( outputPFOMomentum );
//				outputPFO->setEnergy( pfoE );
//				outputPFO->setMass( outputPFO->getMass() );
//				outputPFO->setCharge(outputPFO->getCharge());
				outputPFO->setCovMatrix( outputCovMatrix );
/*				outputPFO->setReferencePoint(outputPFO->getReferencePoint());
				for (unsigned int j=0; j<outputPFO->getParticleIDs().size(); ++j)
				{
					ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>(outputPFO->getParticleIDs()[j]);
				        ParticleIDImpl* outPID = new ParticleIDImpl;
				        outPID->setType(inPID->getType());
				        outPID->setPDG(inPID->getPDG());
				        outPID->setLikelihood(inPID->getLikelihood());
				        outPID->setAlgorithmType(inPID->getAlgorithmType()) ;
				        for (unsigned int k=0; k<inPID->getParameters().size()  ; ++k) outPID->addParameter(inPID->getParameters()[k]) ;
				        outputPFO->addParticleID(outPID);
				}
				outputPFO->setParticleIDUsed(outputPFO->getParticleIDUsed());
				outputPFO->setGoodnessOfPID(outputPFO->getGoodnessOfPID());
				for (unsigned int j=0; j< outputPFO->getParticles().size(); ++j)
				{
					outputPFO->addParticle(outputPFO->getParticles()[j]);
				}
				for (unsigned int j=0; j<outputPFO->getClusters().size(); ++j)
				{
					outputPFO->addCluster(outputPFO->getClusters()[j]);
				}
				for (unsigned int j=0; j<outputPFO->getTracks().size(); ++j)
				{
					outputPFO->addTrack(outputPFO->getTracks()[j]);
				}
				outputPFO->setStartVertex(outputPFO->getStartVertex());
*/
				streamlog_out(DEBUG6) << "	Updated PFO:" << std::endl;
				streamlog_out(DEBUG5) << *outputPFO << std::endl;
			}
			outputPfoCollection->addElement( outputPFO );
		}
		pLCEvent->addCollection( outputPfoCollection , m_outputPfoCollection );
	}
        catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input collection not found in event " << pLCEvent->getEventNumber() << std::endl;
        }

}

void ClusterCovMat::getNeutralCovMatFirstOrder( ReconstructedParticleImpl* NeutralPFO , std::vector<float> CovMat )
{

	Cluster* cluster = dynamic_cast<Cluster*>( NeutralPFO->getClusters()[ 0 ] );
	CalorimeterHitVec hits = cluster->getCalorimeterHits();
	int nhits = cluster->getCalorimeterHits().size();
	std::vector<std::vector<double>> ClusterHitsParameters; ClusterHitsParameters.clear();
	std::vector<double> HitParameters;
	double eCluster = 0.0;
	streamlog_out(DEBUG2) << "Cluster has " << nhits << " hits" << std::endl;
//	streamlog_out(DEBUG2) << 	"n 		E 		x 		y 		z 		" << std::endl;
	for ( int ihit = 0 ; ihit < nhits ; ++ihit )
	{
		CalorimeterHit* calo_hit  = dynamic_cast<CalorimeterHit*>( hits[ ihit ] );
		eCluster += calo_hit->getEnergy();
		HitParameters.clear();
		HitParameters.push_back( calo_hit->getPosition()[ 0 ] );
		HitParameters.push_back( calo_hit->getPosition()[ 1 ] );
		HitParameters.push_back( calo_hit->getPosition()[ 2 ] );
		HitParameters.push_back( calo_hit->getEnergy() );
		ClusterHitsParameters.push_back( HitParameters );
//		streamlog_out(DEBUG2) << 	ihit << "		" << calo_hit->getEnergy() << "	" << calo_hit->getPosition()[ 0 ] << "	" << calo_hit->getPosition()[ 1 ] << "	" << calo_hit->getPosition()[ 2 ] << std::endl;
	}

	double clusterEnergyWeightedMeans[ 4 ]{ 0.0 , 0.0 , 0.0 , 0.0 };
	double weight = 1.0;
	for ( int i_space = 0 ; i_space < 4 ; ++i_space )
	{
		for ( int ihit = 0 ; ihit < nhits ; ++ihit )
		{
			if ( i_space < 3 )
			{
				weight = ClusterHitsParameters[ ihit ][ 3 ] / eCluster;
			}
			else
			{
				weight = 1.0;
			}
			clusterEnergyWeightedMeans[ i_space ] += ClusterHitsParameters[ ihit ][ i_space ] * weight;
		}
	}
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "Cluster Properties" << std::endl;
	streamlog_out(DEBUG2) << 	"E 			x(COG) 		y(COG) 		z(COG) 		" << std::endl;
	streamlog_out(DEBUG2) << 	clusterEnergyWeightedMeans[ 3 ] << "(" << eCluster << ")" << "	" << clusterEnergyWeightedMeans[ 0 ] << "	" << clusterEnergyWeightedMeans[ 1 ] << "	" << clusterEnergyWeightedMeans[ 2 ] << std::endl;

//	Obtain covariance matrix on (px,py,pz,E) from the
//	covariance matrix on cluster parameters (px,py,pz,|p|=Ec).
//	=> E^2 = Ec^2 + m^2	;	|p| = Ec
//	define the jacobian as the 4x4 matrix:
//
//
//
//	Dpx/Dx	= |P|.(r2-x2)/r3	Dpy/Dx	= -|P|.x.y/r3		Dpz/Dx	= -|P|.x.z/r3		DE/Dx	= 0
//
//	Dpx/Dy	= -|P|.y.x/r3		Dpy/Dy	= |P|.(r2-y2)/r3	Dpz/Dy	= -|P|.y.z/r3		DE/Dy	= 0
//
//	Dpx/Dz	= -|P|.z.x/r3		Dpy/Dz	= -|P|.z.y/r3		Dpz/Dz	= |P|.(r2-z2)/r3	DE/Dz	= 0
//
//	Dpx/DEc	= (E/|p|).(x/r)		Dpy/DEc	= (E/|p|).(y/r)		Dpz/DEc	= (E/|p|).(z/r)		DE/DEc	= 1
//


	const int nCluPar = 4;
	const int nPFOPar = 4;
	double X	= clusterEnergyWeightedMeans[ 0 ];
	double Y	= clusterEnergyWeightedMeans[ 1 ];
	double Z	= clusterEnergyWeightedMeans[ 2 ];
	double X2	= pow( clusterEnergyWeightedMeans[ 0 ] , 2 );
	double Y2	= pow( clusterEnergyWeightedMeans[ 1 ] , 2 );
	double Z2	= pow( clusterEnergyWeightedMeans[ 2 ] , 2 );
	double R	= sqrt( pow( clusterEnergyWeightedMeans[ 0 ] , 2 ) + pow( clusterEnergyWeightedMeans[ 1 ] , 2 ) + pow( clusterEnergyWeightedMeans[ 2 ] , 2 ) );
	double R2	= pow( clusterEnergyWeightedMeans[ 0 ] , 2 ) + pow( clusterEnergyWeightedMeans[ 1 ] , 2 ) + pow( clusterEnergyWeightedMeans[ 2 ] , 2 );
	double R3	= R * R2;
	double Momentum = clusterEnergyWeightedMeans[ 3 ];
	double Energy	= sqrt( pow( Momentum , 2 ) + pow( NeutralPFO->getMass() , 2 ) );
	double Px	= Momentum * clusterEnergyWeightedMeans[ 0 ] / R;
	double Py	= Momentum * clusterEnergyWeightedMeans[ 1 ] / R;
	double Pz	= Momentum * clusterEnergyWeightedMeans[ 2 ] / R;
	double D1_PFOPar_CluPar[ nPFOPar ][ nCluPar ];
	D1_PFOPar_CluPar[ 0 ][ 0 ] = Momentum * ( R2 - X2 ) / R3;
	D1_PFOPar_CluPar[ 0 ][ 1 ] = -1.0 * Momentum * X * Y / R3;
	D1_PFOPar_CluPar[ 0 ][ 2 ] = -1.0 * Momentum * X * Z / R3;
	D1_PFOPar_CluPar[ 0 ][ 3 ] = Energy * X / ( Momentum * R );
	D1_PFOPar_CluPar[ 1 ][ 0 ] = -1.0 * Momentum * Y * X / R3;
	D1_PFOPar_CluPar[ 1 ][ 1 ] = Momentum * ( R2 - Y2 ) / R3;
	D1_PFOPar_CluPar[ 1 ][ 2 ] = -1.0 * Momentum * Y * Z / R3;
	D1_PFOPar_CluPar[ 1 ][ 3 ] = Energy * Y / ( Momentum * R );
	D1_PFOPar_CluPar[ 2 ][ 0 ] = -1.0 * Momentum * Z * X / R3;
	D1_PFOPar_CluPar[ 2 ][ 1 ] = -1.0 * Momentum * Z * Y / R3;
	D1_PFOPar_CluPar[ 2 ][ 2 ] = Momentum * ( R2 - Z2 ) / R3;
	D1_PFOPar_CluPar[ 2 ][ 3 ] = Energy * Z / ( Momentum * R );
	D1_PFOPar_CluPar[ 3 ][ 0 ] = 0.0;
	D1_PFOPar_CluPar[ 3 ][ 1 ] = 0.0;
	D1_PFOPar_CluPar[ 3 ][ 2 ] = 0.0;
	D1_PFOPar_CluPar[ 3 ][ 3 ] = 1.0;

	const int rows			= 4; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	TMatrixD Gradient( kspace_time_dim , kspace_time_dim );
	TMatrixD hitMoment( kspace_time_dim , 1 );
	TMatrixD hitMomentTranspose( 1 , kspace_time_dim );
	TMatrixD GradientMoment( kspace_time_dim , 1 );
	TMatrixD MomentTransposeGradientMoment( 1 , 1 );
	double SigmaEhit2 = 0.0;
	for ( int ihit = 0 ; ihit < nhits ; ++ihit ) SigmaEhit2 += pow( ClusterHitsParameters[ ihit ][ 3 ] , 2 );

	double coefficient_FirstOrder = SigmaEhit2 / ( pow( clusterEnergyWeightedMeans[ 3 ] , 2 ) - SigmaEhit2 );
//	double coefficient_FirstOrder = pow( clusterEnergyWeightedMeans[ 3 ] , 2 ) / ( pow( clusterEnergyWeightedMeans[ 3 ] , 2 ) - SigmaEhit2 );
//	double coefficient_FirstOrder = SigmaEhit2 * pow( clusterEnergyWeightedMeans[ 3 ] , 1 ) / ( pow( clusterEnergyWeightedMeans[ 3 ] , 2 ) - SigmaEhit2 ) / pow( clusterEnergyWeightedMeans[ 3 ] , 2 );
//	double coefficient_FirstOrder = 1.0;
	streamlog_out(DEBUG0) << "		sigma(Ehit) = " << clusterEnergyWeightedMeans[ 3 ] << std::endl;
	streamlog_out(DEBUG0) << "		sigma(Ehit^2) = " << SigmaEhit2 << std::endl;
	streamlog_out(DEBUG0) << "		[sigma(Ehit)]^2 = " << pow( clusterEnergyWeightedMeans[ 3 ] , 2 ) << std::endl;

	for ( int ihit = 0 ; ihit < nhits ; ++ihit )
	{
		for ( int index = 0 ; index < 4 ; ++index )
		{
			hitMoment( index , 0 ) = ClusterHitsParameters[ ihit ][ index ] - clusterEnergyWeightedMeans[ index ];
			hitMomentTranspose( 0 , index ) = ClusterHitsParameters[ ihit ][ index ] - clusterEnergyWeightedMeans[ index ];
		}
		for ( int i = 0 ; i < kspace_time_dim ; ++i )
		{
			for ( int j = 0 ; j < kspace_time_dim ; ++j )
			{
				TMatrixD gradientPar1( kspace_time_dim , 1 );
				TMatrixD gradientPar2( 1 , kspace_time_dim );
				for ( int index = 0 ; index < kspace_time_dim ; ++index )
				{
					gradientPar1( index , 0 ) = D1_PFOPar_CluPar[ i ][ index ];
					gradientPar2( 0 , index ) = D1_PFOPar_CluPar[ j ][ index ];
				}
				Gradient.Mult( gradientPar1 , gradientPar2 );
				GradientMoment.Mult( Gradient , hitMoment );
				MomentTransposeGradientMoment.Mult( hitMomentTranspose , GradientMoment );
				covMatrixMomenta( i , j ) += coefficient_FirstOrder * MomentTransposeGradientMoment( 0 , 0 ) * ClusterHitsParameters[ ihit ][ 3 ] / clusterEnergyWeightedMeans[ 3 ];
			}
		}
	}

	double dPhi_dPx		= -Py / ( pow( Px , 2 ) + pow( Py , 2 ) );
	double dPhi_dPy		= Px / ( pow( Px , 2 ) + pow( Py , 2 ) );
	double dPhi_dx		= -Y / ( pow( X , 2 ) + pow( Y , 2 ) );
	double dPhi_dy		= X / ( pow( X , 2 ) + pow( Y , 2 ) );



	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG0) << "		Input PFO Covariance Matrix:" << std::endl;
	streamlog_out(DEBUG0) << "			" << NeutralPFO->getCovMatrix()[ 0 ] << "	, " << NeutralPFO->getCovMatrix()[ 1 ] << "	, " << NeutralPFO->getCovMatrix()[ 3 ] << "	, " << NeutralPFO->getCovMatrix()[ 6 ] << std::endl;
	streamlog_out(DEBUG0) << "			" << NeutralPFO->getCovMatrix()[ 1 ] << "	, " << NeutralPFO->getCovMatrix()[ 2 ] << "	, " << NeutralPFO->getCovMatrix()[ 4 ] << "	, " << NeutralPFO->getCovMatrix()[ 7 ] << std::endl;
	streamlog_out(DEBUG0) << "			" << NeutralPFO->getCovMatrix()[ 3 ] << "	, " << NeutralPFO->getCovMatrix()[ 4 ] << "	, " << NeutralPFO->getCovMatrix()[ 5 ] << "	, " << NeutralPFO->getCovMatrix()[ 8 ] << std::endl;
	streamlog_out(DEBUG0) << "			" << NeutralPFO->getCovMatrix()[ 6 ] << "	, " << NeutralPFO->getCovMatrix()[ 7 ] << "	, " << NeutralPFO->getCovMatrix()[ 8 ] << "	, " << NeutralPFO->getCovMatrix()[ 9 ] << std::endl;
	streamlog_out(DEBUG0) << "	sigmaPhi2 = " << pow( dPhi_dPx , 2 ) * NeutralPFO->getCovMatrix()[ 0 ] + pow( dPhi_dPy , 2 ) * NeutralPFO->getCovMatrix()[ 2 ] + 2.0 * dPhi_dPx * dPhi_dPy * NeutralPFO->getCovMatrix()[ 1 ] << std::endl;

	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG0) << "		Output PFO Covariance Matrix:" << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrixMomenta( 0 , 0 ) << "	, " << covMatrixMomenta( 0 , 1 ) << "	, " << covMatrixMomenta( 0 , 2 ) << "	, " << covMatrixMomenta( 0 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrixMomenta( 1 , 0 ) << "	, " << covMatrixMomenta( 1 , 1 ) << "	, " << covMatrixMomenta( 1 , 2 ) << "	, " << covMatrixMomenta( 1 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrixMomenta( 2 , 0 ) << "	, " << covMatrixMomenta( 2 , 1 ) << "	, " << covMatrixMomenta( 2 , 2 ) << "	, " << covMatrixMomenta( 2 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrixMomenta( 3 , 0 ) << "	, " << covMatrixMomenta( 3 , 1 ) << "	, " << covMatrixMomenta( 3 , 2 ) << "	, " << covMatrixMomenta( 3 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "	sigmaPhi2 = " << pow( dPhi_dPx , 2 ) * covMatrixMomenta( 0 , 0 ) + pow( dPhi_dPy , 2 ) * covMatrixMomenta( 1 , 1 ) + 2.0 * dPhi_dPx * dPhi_dPy * covMatrixMomenta( 1 , 0 ) << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG0) << "		Ratio" << std::endl;
	streamlog_out(DEBUG0) << "			" << NeutralPFO->getCovMatrix()[ 0 ] / covMatrixMomenta( 0 , 0 ) << "	, " << NeutralPFO->getCovMatrix()[ 1 ] / covMatrixMomenta( 0 , 1 ) << "	, " << NeutralPFO->getCovMatrix()[ 3 ] / covMatrixMomenta( 0 , 2 ) << "	, " << NeutralPFO->getCovMatrix()[ 6 ] / covMatrixMomenta( 0 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << NeutralPFO->getCovMatrix()[ 1 ] / covMatrixMomenta( 1 , 0 ) << "	, " << NeutralPFO->getCovMatrix()[ 2 ] / covMatrixMomenta( 1 , 1 ) << "	, " << NeutralPFO->getCovMatrix()[ 4 ] / covMatrixMomenta( 1 , 2 ) << "	, " << NeutralPFO->getCovMatrix()[ 7 ] / covMatrixMomenta( 1 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << NeutralPFO->getCovMatrix()[ 3 ] / covMatrixMomenta( 2 , 0 ) << "	, " << NeutralPFO->getCovMatrix()[ 4 ] / covMatrixMomenta( 2 , 1 ) << "	, " << NeutralPFO->getCovMatrix()[ 5 ] / covMatrixMomenta( 2 , 2 ) << "	, " << NeutralPFO->getCovMatrix()[ 8 ] / covMatrixMomenta( 2 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << NeutralPFO->getCovMatrix()[ 6 ] / covMatrixMomenta( 3 , 0 ) << "	, " << NeutralPFO->getCovMatrix()[ 7 ] / covMatrixMomenta( 3 , 1 ) << "	, " << NeutralPFO->getCovMatrix()[ 8 ] / covMatrixMomenta( 3 , 2 ) << "	, " << NeutralPFO->getCovMatrix()[ 9 ] / covMatrixMomenta( 3 , 3 ) << std::endl;
	CovMat.clear();



//	CovMat.push_back( covMatrixMomenta(0,0) ); // x-x
//	CovMat.push_back( covMatrixMomenta(1,0) ); // y-x
//	CovMat.push_back( covMatrixMomenta(1,1) ); // y-y
//	CovMat.push_back( covMatrixMomenta(2,0) ); // z-x
//	CovMat.push_back( covMatrixMomenta(2,1) ); // z-y
//	CovMat.push_back( covMatrixMomenta(2,2) ); // z-z
//	CovMat.push_back( covMatrixMomenta(3,0) ); // e-x
//	CovMat.push_back( covMatrixMomenta(3,1) ); // e-y
//	CovMat.push_back( covMatrixMomenta(3,2) ); // e-z
//	CovMat.push_back( covMatrixMomenta(3,3) ); // e-e
	CovMat.push_back( 0.0 ); // x-x
	CovMat.push_back( 0.0 ); // y-x
	CovMat.push_back( 0.0 ); // y-y
	CovMat.push_back( 0.0 ); // z-x
	CovMat.push_back( 0.0 ); // z-y
	CovMat.push_back( 0.0 ); // z-z
	CovMat.push_back( 0.0 ); // e-x
	CovMat.push_back( 0.0 ); // e-y
	CovMat.push_back( 0.0 ); // e-z
	CovMat.push_back( 0.0 ); // e-e
	if ( CovMat.size() == 10 ) streamlog_out(DEBUG0) << "	FourMomentumCovarianceMatrix Filled succesfully" << std::endl;

//	return CovMat;

}

void ClusterCovMat::getNeutralCovMatSecondOrder( ReconstructedParticleImpl* NeutralPFO , std::vector<float> CovMat )
{

	Cluster* cluster = dynamic_cast<Cluster*>( NeutralPFO->getClusters()[ 0 ] );
	CalorimeterHitVec hits = cluster->getCalorimeterHits();
	int nhits = cluster->getCalorimeterHits().size();
	std::vector<double> ehit,xhit,yhit,zhit;
	for ( int ihit = 0 ; ihit < nhits ; ++ihit )
	{
		CalorimeterHit* calo_hit  = dynamic_cast<CalorimeterHit*>( hits[ ihit ] );
		ehit.push_back( calo_hit->getEnergy() );
		xhit.push_back( calo_hit->getPosition()[ 0 ] );
		yhit.push_back( calo_hit->getPosition()[ 1 ] );
		zhit.push_back( calo_hit->getPosition()[ 2 ] );
	}

	const int rows			= 4; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP; covP.clear();


}

std::vector<float> ClusterCovMat::getNeutralCovMatFirstOrder( TVector3 clusterPosition , float pfoEc , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError )
{

//	Obtain covariance matrix on (px,py,pz,E) from the
//	covariance matrix on cluster parameters (px,py,pz,|p|=Ec).
//	=> E^2 = Ec^2 + m^2	;	|p| = Ec
//	define the jacobian as the 4x4 matrix:
//
//
//
//			Dpx/Dx			Dpy/Dx			Dpz/Dx			DE/Dx
//
//			Dpx/Dy			Dpy/Dy			Dpz/Dy			DE/Dy
//	J =
//			Dpx/Dz			Dpy/Dz			Dpz/Dz			DE/Dz
//
//			Dpx/DEc			Dpy/DEc			Dpz/DEc			DE/DEc
//
//
//
//
//
//			 |P|.(r2-x2)/r3		-|P|.x.y/r3		-|P|.x.z/r3		0
//
//			-|P|.y.x/r3		 |P|.(r2-y2)/r3		-|P|.y.z/r3		0
//	J =
//			-|P|.z.x/r3		-|P|.z.y/r3		 |P|.(r2-z2)/r3		0
//
//			 (E/|p|).(x/r)		 (E/|p|).(y/r)		 (E/|p|).(z/r)		1
//
//
//
//
//	CovMatrix elements in terms of cluster position error and cluster energy error:
//
//			x.x			x.y			x.z			x.Ec
//
//			y.x			y.y			y.z			y.Ec
//	Cov =
//			z.x			z.y			z.z			z.Ec
//
//			Ec.x			Ec.y			Ec.z			Ec.Ec
//
//
//

	const int rows			= 4; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP; covP.clear();

//	pfoMass			= 0.0;

	float pfoX		=	clusterPosition.X();
	float pfoY		=	clusterPosition.Y();
	float pfoZ		=	clusterPosition.Z();
	float pfoR		=	std::sqrt( pow( pfoX , 2 ) + pow( pfoY , 2 ) + pow( pfoZ , 2 ) );
	float pfoX2		=	pow( pfoX , 2 );
	float pfoY2		=	pow( pfoY , 2 );
	float pfoZ2		=	pow( pfoZ , 2 );
	float pfoR2		=	pow( pfoR , 2 );
	float pfoR3		=	pow( pfoR , 3 );
	float SigmaX2		=	clusterPositionError[ 0 ];
	float SigmaXY		=	clusterPositionError[ 1 ];
	float SigmaY2		=	clusterPositionError[ 2 ];
	float SigmaXZ		=	clusterPositionError[ 3 ];
	float SigmaYZ		=	clusterPositionError[ 4 ];
	float SigmaZ2		=	clusterPositionError[ 5 ];
	float SigmaE2		=	pow( clusterEnergyError , 2 );

	float pfoP = ( m_isClusterEnergyKinEnergy ? sqrt( pow( pfoEc , 2 ) + 2 * pfoMass * pfoEc ) : pfoEc );
	float pfoE = ( m_isClusterEnergyKinEnergy ? sqrt( pow( pfoP , 2 ) + pow( pfoMass , 2 ) ) : sqrt( pow( pfoP , 2 ) + pow( pfoMass , 2 ) ) );
	float derivative_coeff	= ( ( m_useTrueJacobian && !m_isClusterEnergyKinEnergy ) ? pfoP / pfoE : 1.0 );

	streamlog_out(DEBUG0) << "	Cluster information obtained:" << std::endl;
	streamlog_out(DEBUG0) << "		Cluster Information:" << std::endl;
	streamlog_out(DEBUG0) << "		X = " << pfoX << "	, Y = " << pfoY << "	, Z = " << pfoZ << "	, Energy = " << pfoEc << "	, Mass = " << pfoMass << std::endl;
//	cluster covariance matrix by rows
	double cluster_cov_matrix_by_rows[rows*rows] =
	{
		SigmaX2		,	SigmaXY		,	SigmaXZ		,	0	,
		SigmaXY		,	SigmaY2		,	SigmaYZ		,	0	,
		SigmaXZ		,	SigmaYZ		,	SigmaZ2		,	0	,
		0		,	0		,	0		,	SigmaE2
	};
	TMatrixD covMatrix_cluster(rows,rows, cluster_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG0) << "	Cluster covariance matrix array converted to cluster covariance matrix" << std::endl;

	streamlog_out(DEBUG0) << "		Cluster Position Error:" << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrix_cluster( 0 , 0 ) << "	, " << covMatrix_cluster( 0 , 1 ) << "	, " << covMatrix_cluster( 0 , 2 ) << "	, " << covMatrix_cluster( 0 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrix_cluster( 1 , 0 ) << "	, " << covMatrix_cluster( 1 , 1 ) << "	, " << covMatrix_cluster( 1 , 2 ) << "	, " << covMatrix_cluster( 1 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrix_cluster( 2 , 0 ) << "	, " << covMatrix_cluster( 2 , 1 ) << "	, " << covMatrix_cluster( 2 , 2 ) << "	, " << covMatrix_cluster( 2 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrix_cluster( 3 , 0 ) << "	, " << covMatrix_cluster( 3 , 1 ) << "	, " << covMatrix_cluster( 3 , 2 ) << "	, " << covMatrix_cluster( 3 , 3 ) << std::endl;


//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		pfoP * ( pfoR2 - pfoX2 ) / pfoR3			,	-pfoP * pfoX * pfoY / pfoR3				,	-pfoP * pfoX * pfoZ / pfoR3				,	0			,
		-pfoP * pfoY * pfoX / pfoR3				,	pfoP * ( pfoR2 - pfoY2 ) / pfoR3			,	-pfoP * pfoY * pfoZ / pfoR3				,	0			,
		-pfoP * pfoZ * pfoX / pfoR3				,	-pfoP * pfoZ * pfoY / pfoR3				,	pfoP * ( pfoR2 - pfoZ2 ) / pfoR3			,	0			,
		derivative_coeff * pfoE * pfoX / ( pfoP * pfoR )	,	derivative_coeff * pfoE * pfoY / ( pfoP * pfoR )	,	derivative_coeff * pfoE * pfoZ / ( pfoP * pfoR )	,	derivative_coeff
	};

//	construct the Jacobian using previous array ("F" if filling by columns, "C" if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG0) << "	Jacobian array converted to Jacobian matrix" << std::endl;
	streamlog_out(DEBUG0) << "		Jacobian:" << std::endl;
	streamlog_out(DEBUG0) << "			" << jacobian( 0 , 0 ) << "	, " << jacobian( 0 , 1 ) << "	, " << jacobian( 0 , 2 ) << "	, " << jacobian( 0 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << jacobian( 1 , 0 ) << "	, " << jacobian( 1 , 1 ) << "	, " << jacobian( 1 , 2 ) << "	, " << jacobian( 1 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << jacobian( 2 , 0 ) << "	, " << jacobian( 2 , 1 ) << "	, " << jacobian( 2 , 2 ) << "	, " << jacobian( 2 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << jacobian( 3 , 0 ) << "	, " << jacobian( 3 , 1 ) << "	, " << jacobian( 3 , 2 ) << "	, " << jacobian( 3 , 3 ) << std::endl;


	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_cluster) ,
					jacobian
					);

	streamlog_out(DEBUG0) << "		PFO Covariance Matrix:" << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrixMomenta( 0 , 0 ) << "	, " << covMatrixMomenta( 0 , 1 ) << "	, " << covMatrixMomenta( 0 , 2 ) << "	, " << covMatrixMomenta( 0 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrixMomenta( 1 , 0 ) << "	, " << covMatrixMomenta( 1 , 1 ) << "	, " << covMatrixMomenta( 1 , 2 ) << "	, " << covMatrixMomenta( 1 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrixMomenta( 2 , 0 ) << "	, " << covMatrixMomenta( 2 , 1 ) << "	, " << covMatrixMomenta( 2 , 2 ) << "	, " << covMatrixMomenta( 2 , 3 ) << std::endl;
	streamlog_out(DEBUG0) << "			" << covMatrixMomenta( 3 , 0 ) << "	, " << covMatrixMomenta( 3 , 1 ) << "	, " << covMatrixMomenta( 3 , 2 ) << "	, " << covMatrixMomenta( 3 , 3 ) << std::endl;

	covP.push_back( covMatrixMomenta(0,0) ); // x-x
	covP.push_back( covMatrixMomenta(1,0) ); // y-x
	covP.push_back( covMatrixMomenta(1,1) ); // y-y
	covP.push_back( covMatrixMomenta(2,0) ); // z-x
	covP.push_back( covMatrixMomenta(2,1) ); // z-y
	covP.push_back( covMatrixMomenta(2,2) ); // z-z
	covP.push_back( covMatrixMomenta(3,0) ); // e-x
	covP.push_back( covMatrixMomenta(3,1) ); // e-y
	covP.push_back( covMatrixMomenta(3,2) ); // e-z
	covP.push_back( covMatrixMomenta(3,3) ); // e-e
	if ( covP.size() == 10 ) streamlog_out(DEBUG0) << "	FourMomentumCovarianceMatrix Filled succesfully" << std::endl;

	return covP;

}

void ClusterCovMat::check(EVENT::LCEvent *pLCEvent)
{

	LCCollection *inputPfoCollection{};
	LCCollection *outputPfoCollection{};
	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		outputPfoCollection = pLCEvent->getCollection(m_outputPfoCollection);
		int n_inputPFOs = inputPfoCollection->getNumberOfElements();
		int n_outputPFOs = outputPfoCollection->getNumberOfElements();
		streamlog_out(DEBUG) << " CHECK : processed event: " << pLCEvent->getEventNumber() << " (Number of inputPFOS: " << n_inputPFOs << " , Number of outputPFOs: " << n_outputPFOs <<")" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input/Output collection not found in event " << pLCEvent->getEventNumber() << std::endl;
        }

}

void ClusterCovMat::end()
{

//	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}
