package dynamicprogramming;

import base.Array;
import model.Round0Model;
import base.Sequence;
import base.Shape;

public class FullFeatureNucleotideDinucleotideShape implements DynamicProgramming {
	private boolean isNSBinding;
	private int l;
	private int k;
	private int r0k;
	private int flankLength;
	private int maxFrames;
	private int maxSubFeatures;
	private int nShapeClasses;
	private int totNucFeatures;
	private int totDinucFeatures;
	private int totShapeFeatures;
	private int totFeatures;
	private long leftFlank;
	private long leftFlankRC;
	private long rightFlank;
	private long rightFlankRC;
	private long subFeatureMask;
	private long r0FeatureMask;
	private double Z;
	private double[] nucAlphas;
	private double[] nucGradients;
	private double[] dinucAlphas;
	private double[] dinucGradients;
	private double[] shapeBetas;
	private double[] shapeGradients;
	private double[] r0Alphas;
	private double[][] shapeFeatures;
	private double[][] fsInitMatrix;
	private double[][] rsInitMatrix;
	private double[][] hessian;
	
	public FullFeatureNucleotideDinucleotideShape(int l, int k, boolean isNSBinding, int flankLength, 
			String lFlank, String rFlank, Round0Model R0Model, Shape shapeModel) {
		this.l 				= l;
		this.k 				= k;
		this.isNSBinding	= isNSBinding;
		this.flankLength	= flankLength;
		r0k 				= R0Model.getK();
		maxFrames 			= l-k+1+2*flankLength;
		maxSubFeatures		= (int) Math.max(256, Math.pow(4, r0k-1));
		shapeFeatures 		= shapeModel.getFeatures();
		nShapeClasses		= shapeModel.nShapeFeatures();
		totNucFeatures		= 4*k;
		totDinucFeatures	= 16*(k-1);
		totShapeFeatures	= nShapeClasses*k;
		totFeatures			= totNucFeatures + totDinucFeatures + totShapeFeatures;
		r0Alphas			= R0Model.getAlphas();
		leftFlank 			= (new Sequence(lFlank, 0, lFlank.length())).getValue();
		rightFlankRC 		= reverseComplement(leftFlank, lFlank.length());
		rightFlankRC 		= reverse(rightFlankRC, lFlank.length());
		rightFlank 			= (new Sequence(rFlank, 0, rFlank.length())).getValue();
		leftFlankRC 		= reverseComplement(rightFlank, rFlank.length());
		rightFlank 			= reverse(rightFlank, rFlank.length());
		subFeatureMask		= (long) maxSubFeatures - 1;
		r0FeatureMask		= (long) Math.pow(4, r0k) - 1;				//Mask to select the Round0 Feature
		long flankMask		= (long) Math.pow(4, flankLength) - 1;
				
		int fsLFlankInit	= (int) ((leftFlank >> 2*flankLength) & (long) (Math.pow(4, Math.max(r0k-1-flankLength, 2)) - 1));
		int rsLFlankInit	= (int) ((leftFlankRC >> 2*flankLength) & (long) (Math.pow(4, Math.max(r0k-1-flankLength, 2)) - 1));
		
		leftFlank			= leftFlank & flankMask;
		leftFlank			= reverse(leftFlank, flankLength);
		leftFlankRC			= leftFlankRC & flankMask;
		leftFlankRC			= reverse(leftFlankRC, flankLength);
		
		//Create initialization matrix
		fsInitMatrix 		= new double[maxFrames][maxSubFeatures];
		rsInitMatrix		= new double[maxFrames][maxSubFeatures];
		for (int i=0; i<maxFrames; i++) {
			fsInitMatrix[i][(int) fsLFlankInit] = 1;
			rsInitMatrix[i][(int) rsLFlankInit] = 1;
		}		
	}
	
	@Override
	public void setAlphas(double[] nucAlphas, double[] dinucAlphas, double[] shapeBetas) {
		this.nucAlphas 		= nucAlphas;
		this.dinucAlphas 	= dinucAlphas;
		this.shapeBetas		= shapeBetas;
	}
	
	public double recursiveZ() {
		int startOffset;
		int endOffset;
		int shapeStartOffset;
		int shapeEndOffset;
		int position;
		int shapePosition;
		int fShapeIndex;
		int rShapeIndex;
		double fR0Weight	  	= 0;
		double rR0Weight		= 0;
		double fNucWeight		= 0;
		double rNucWeight		= 0;
		double fDinucWeight	  	= 0;
		double rDinucWeight		= 0;
		double fShapeWeight		= 0;
		double rShapeWeight		= 0;
		double[][] fsOffsets 	= matrixClone(fsInitMatrix); 			//Matrix to store current and previous dynamic programming values
		double[][] rsOffsets	= matrixClone(rsInitMatrix);
		double[][] fsWork		= new double[maxFrames][maxSubFeatures];
		double[][] rsWork		= new double[maxFrames][maxSubFeatures];
		long fNewWorkFeat;
		long rNewWorkFeat;
		long fsLFlank			= leftFlank;
		long fsRFlank			= rightFlank;
		long rsLFlank			= leftFlankRC;
		long rsRFlank			= rightFlankRC;
		
		//Loop over left fixed region
		for (int index=0; index<flankLength; index++) {									//NO ROUND0 MODEL HERE
			startOffset 	 = (index-k+1 < 0) ? 0 : index-k+1;
			endOffset 		 = (index > maxFrames-1) ? maxFrames-1 : index;
			shapeStartOffset = (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset 	 = (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			
			//Find new base
			int fNewBase = (int) fsLFlank & 3;
			int rNewBase = (int) rsLFlank & 3;
			fsLFlank >>= 2;
			rsLFlank >>= 2;
			for (long currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				fNewWorkFeat = (currSubFeat<<2) | fNewBase;
				rNewWorkFeat = (currSubFeat<<2) | rNewBase;
				for (int currOffset=0; currOffset<maxFrames; currOffset++) {
					position 		= index-currOffset;
					shapePosition 	= position-2;
					if (currOffset>=startOffset && currOffset<=endOffset) {
						if (position < 1) {
							fDinucWeight = 1;
							rDinucWeight = 1;
						} else {
							fDinucWeight = dinucAlphas[(position-1)*16 + (int) (fNewWorkFeat & 15)];
							rDinucWeight = dinucAlphas[(position-1)*16 + (int) (rNewWorkFeat & 15)];
						}
						fNucWeight = nucAlphas[position*4 + fNewBase];
						rNucWeight = nucAlphas[position*4 + rNewBase];
					} else {
						fDinucWeight = 1;
						rDinucWeight = 1;
						fNucWeight = 1;
						rNucWeight = 1;
					}
					if (index >= 2 && currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						fShapeIndex	= (int) (fNewWorkFeat & 1023);
						rShapeIndex = (int) (rNewWorkFeat & 1023);
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosition*nShapeClasses + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosition*nShapeClasses + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fShapeWeight = Math.exp(fShapeWeight);
						rShapeWeight = Math.exp(rShapeWeight);
					} else {
						fShapeWeight = 1;
						rShapeWeight = 1;
					}
					fsWork[currOffset][(int) (fNewWorkFeat & subFeatureMask)] += 
							fsOffsets[currOffset][(int) currSubFeat]*fNucWeight*fDinucWeight*fShapeWeight;
					rsWork[currOffset][(int) (rNewWorkFeat & subFeatureMask)] += 
							rsOffsets[currOffset][(int) currSubFeat]*rNucWeight*rDinucWeight*rShapeWeight;
				}
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
		}
		
		//Loop over variable region
		for (int index=flankLength; index<l+flankLength; index++) {								//R0, NUC, DINUC, SHAPE
			startOffset 	 = (index-k+1 < 0) ? 0 : index-k+1;
			endOffset 		 = (index > maxFrames-1) ? maxFrames-1 : index;
			shapeStartOffset = (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset 	 = (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			for (long currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				for (int newBase=0; newBase<4; newBase++) {
					fNewWorkFeat = (currSubFeat<<2) | newBase;
					rNewWorkFeat = (currSubFeat<<2) | newBase;
					fR0Weight = r0Alphas[(int) (fNewWorkFeat & r0FeatureMask)];
					rR0Weight = r0Alphas[(int) reverseComplement((rNewWorkFeat & r0FeatureMask), r0k)];
					for (int currOffset=0; currOffset<maxFrames; currOffset++) {
						position = index-currOffset;
						shapePosition = position-2;
						if (currOffset>=startOffset && currOffset<=endOffset) {
							if (position < 1) {
								fDinucWeight = 1;
								rDinucWeight = 1;
							} else {
								fDinucWeight = dinucAlphas[(position-1)*16 + (int) (fNewWorkFeat & 15)];
								rDinucWeight = dinucAlphas[(position-1)*16 + (int) (rNewWorkFeat & 15)];								
							}
							fNucWeight = nucAlphas[position*4 + newBase];
							rNucWeight = nucAlphas[position*4 + newBase];
						} else {
							fDinucWeight = 1;
							rDinucWeight = 1;
							fNucWeight = 1;
							rNucWeight = 1;
						}
						if (index >= 2 && currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							fShapeWeight = 0;
							rShapeWeight = 0;
							fShapeIndex = (int) (fNewWorkFeat & 1023);
							rShapeIndex = (int) (rNewWorkFeat & 1023);
							for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
								fShapeWeight += shapeBetas[shapePosition*nShapeClasses + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
								rShapeWeight += shapeBetas[shapePosition*nShapeClasses + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
							}
							fShapeWeight = Math.exp(fShapeWeight);
							rShapeWeight = Math.exp(rShapeWeight);
						} else {
							fShapeWeight = 1;
							rShapeWeight = 1;
						}
						fsWork[currOffset][(int) (fNewWorkFeat & subFeatureMask)] += 
								fsOffsets[currOffset][(int) currSubFeat]*fR0Weight*fNucWeight*fDinucWeight*fShapeWeight;
						rsWork[currOffset][(int) (rNewWorkFeat & subFeatureMask)] += 
								rsOffsets[currOffset][(int) currSubFeat]*rR0Weight*rNucWeight*rDinucWeight*rShapeWeight;
					}
				}
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
		}
		//Loop over right fixed region
		for (int index=l+flankLength; index<l+2*flankLength; index++) {								//R0, NUC, DINUC, SHAPE
			startOffset 	 = (index-k+1 < 0) ? 0 : index-k+1;
			endOffset 		 = (index > maxFrames-1) ? maxFrames-1 : index;
			shapeStartOffset = (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset 	 = (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			
			//Find new base
			int fNewBase = (int) fsRFlank & 3;
			int rNewBase = (int) rsRFlank & 3;
			fsRFlank >>= 2;
			rsRFlank >>= 2;
			for (long currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				fNewWorkFeat = (currSubFeat<<2) | fNewBase;
				rNewWorkFeat = (currSubFeat<<2) | rNewBase;
				if (index >= l+flankLength+r0k-1) {
					fR0Weight = 1;
					rR0Weight = 1;
				} else {
					fR0Weight = r0Alphas[(int) (fNewWorkFeat & r0FeatureMask)];
					rR0Weight = r0Alphas[(int) reverseComplement((rNewWorkFeat & r0FeatureMask), r0k)];
				}
				for (int currOffset=0; currOffset<maxFrames; currOffset++) {
					position = index-currOffset;
					shapePosition = position-2;
					if (currOffset>=startOffset && currOffset<=endOffset) {
						if (position < 1) {		//dont think the IF statement here is necessary
							fDinucWeight = 1;
							rDinucWeight = 1;
						} else {
							fDinucWeight = dinucAlphas[(position-1)*16 + (int) (fNewWorkFeat & 15)];		
							rDinucWeight = dinucAlphas[(position-1)*16 + (int) (rNewWorkFeat & 15)];			
						}
						fNucWeight = nucAlphas[position*4 + fNewBase];
						rNucWeight = nucAlphas[position*4 + rNewBase];
					} else {
						fDinucWeight = 1;
						rDinucWeight = 1;
						fNucWeight = 1;
						rNucWeight = 1;
					}
					if (index >= 2 && currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						fShapeIndex = (int) (fNewWorkFeat & 1023);
						rShapeIndex = (int) (rNewWorkFeat & 1023);
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosition*nShapeClasses + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosition*nShapeClasses + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fShapeWeight = Math.exp(fShapeWeight);
						rShapeWeight = Math.exp(rShapeWeight);
					} else {
						fShapeWeight = 1;
						rShapeWeight = 1;
					}
					fsWork[currOffset][(int) (fNewWorkFeat & subFeatureMask)] += 
							fsOffsets[currOffset][(int) currSubFeat]*fR0Weight*fNucWeight*fDinucWeight*fShapeWeight;
					rsWork[currOffset][(int) (rNewWorkFeat & subFeatureMask)] += 
							rsOffsets[currOffset][(int) currSubFeat]*rR0Weight*rNucWeight*rDinucWeight*rShapeWeight;
				}
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
		}
		//Loop over last 2 bases to complete shape features
		int maxIndex = Math.max(2, r0k-1-flankLength);
		for (int index=l+2*flankLength; index<l+2*flankLength+maxIndex; index++) {	
			shapeStartOffset = (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset 	 = (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			
			//Find new base
			int fNewBase = (int) fsRFlank & 3;
			int rNewBase = (int) rsRFlank & 3;
			fsRFlank >>= 2;
			rsRFlank >>= 2;
			for (long currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				fNewWorkFeat = (currSubFeat<<2) | fNewBase;			//What is the new feature?
				rNewWorkFeat = (currSubFeat<<2) | rNewBase;
				if (index >= l+flankLength+r0k-1) {
					fR0Weight = 1;
					rR0Weight = 1;
				} else {
					fR0Weight = r0Alphas[(int) (fNewWorkFeat & r0FeatureMask)];
					rR0Weight = r0Alphas[(int) reverseComplement((rNewWorkFeat & r0FeatureMask), r0k)];					
				}
				for (int currOffset=0; currOffset<maxFrames; currOffset++) {
					shapePosition = index-currOffset-2;
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						fShapeIndex = (int) (fNewWorkFeat & 1023);
						rShapeIndex = (int) (rNewWorkFeat & 1023);
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosition*nShapeClasses + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosition*nShapeClasses + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fShapeWeight = Math.exp(fShapeWeight);
						rShapeWeight = Math.exp(rShapeWeight);
					} else {
						fShapeWeight = 1;
						rShapeWeight = 1;
					}
					fsWork[currOffset][(int) (fNewWorkFeat & subFeatureMask)] += 
							fsOffsets[currOffset][(int) currSubFeat]*fShapeWeight*fR0Weight;
					rsWork[currOffset][(int) (rNewWorkFeat & subFeatureMask)] += 
							rsOffsets[currOffset][(int) currSubFeat]*rShapeWeight*rR0Weight;
				}				
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
		}
		
		Z = 0;
		for (int i=0; i<maxFrames; i++) {
			Z += Array.sum(fsOffsets[i]) + Array.sum(rsOffsets[i]);
		}
		
		fsOffsets	= null;
		rsOffsets	= null;
		fsWork		= null;
		rsWork		= null;
		
		return Z;
	}

	@Override
	public void recursiveGradient() {
		int startOffset;
		int endOffset;
		int shapeStartOffset;
		int shapeEndOffset;
		int position;
		int nucPosOffset;
		int fDinucPosOffset;
		int rDinucPosOffset;
		int shapePosOffset;
		int fNextSubFeat;
		int rNextSubFeat;
		int fShapeIndex					= 0;
		int rShapeIndex					= 0;
		double fWeight					= 0;
		double rWeight					= 0;
		double fShapeWeight				= 0;
		double rShapeWeight				= 0;
		double fR0Weight				= 0;
		double rR0Weight				= 0;
		double[][] fsOffsets			= matrixClone(fsInitMatrix);
		double[][] rsOffsets			= matrixClone(rsInitMatrix);
		double[][] fsWork				= new double[maxFrames][maxSubFeatures];
		double[][] rsWork				= new double[maxFrames][maxSubFeatures];
		double[][][] fsNucGradients		= new double[totNucFeatures][maxFrames][maxSubFeatures];
		double[][][] rsNucGradients		= new double[totNucFeatures][maxFrames][maxSubFeatures];
		double[][][] fsDinucGradients	= new double[totDinucFeatures][maxFrames][maxSubFeatures];
		double[][][] rsDinucGradients	= new double[totDinucFeatures][maxFrames][maxSubFeatures];
		double[][][] fsShapeGradients	= new double[totShapeFeatures][maxFrames][maxSubFeatures];
		double[][][] rsShapeGradients	= new double[totShapeFeatures][maxFrames][maxSubFeatures];
		double[][][] fsNucGradientsWork = new double[totNucFeatures][maxFrames][maxSubFeatures];
		double[][][] rsNucGradientsWork = new double[totNucFeatures][maxFrames][maxSubFeatures];
		double[][][] fsDinucGradientsWork=new double[totDinucFeatures][maxFrames][maxSubFeatures];
		double[][][] rsDinucGradientsWork=new double[totDinucFeatures][maxFrames][maxSubFeatures];
		double[][][] fsShapeGradientsWork=new double[totShapeFeatures][maxFrames][maxSubFeatures];
		double[][][] rsShapeGradientsWork=new double[totShapeFeatures][maxFrames][maxSubFeatures];
		long fNewWorkFeat;
		long rNewWorkFeat;
		long fsLFlank					= leftFlank;
		long fsRFlank					= rightFlank;
		long rsLFlank					= leftFlankRC;
		long rsRFlank					= rightFlankRC;

		//Loop over left fixed region
		for (int index=0; index<flankLength; index++) {									//NO ROUND0 MODEL HERE
			startOffset			= (index-k+1 < 0) ? 0 : index-k+1;
			endOffset			= (index > maxFrames-1) ? maxFrames-1 : index;
			shapeStartOffset	= (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset		= (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradientsWork[k][i][j] = 0;
						rsNucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradientsWork[k][i][j] = 0;
						rsDinucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradientsWork[k][i][j] = 0;
						rsShapeGradientsWork[k][i][j] = 0;
					}
				}
			}
			
			//Find new base
			int fNewBase = (int) fsLFlank & 3;
			int rNewBase = (int) rsLFlank & 3;
			fsLFlank >>= 2;
			rsLFlank >>= 2;
			for (int currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				fNewWorkFeat= (currSubFeat<<2) | fNewBase;
				rNewWorkFeat= (currSubFeat<<2) | rNewBase;
				fNextSubFeat= (int) (fNewWorkFeat & subFeatureMask);
				rNextSubFeat= (int) (rNewWorkFeat & subFeatureMask);
				fShapeIndex	= (int) (fNewWorkFeat & 1023);
				rShapeIndex = (int) (rNewWorkFeat & 1023);
				
				for (int currOffset=0; currOffset<startOffset; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset 	= (index-currOffset-2)*nShapeClasses;
					fWeight			= 1;
					rWeight			= 1;
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					//First, update gradient features
					for (int cgf=0; cgf<totNucFeatures; cgf++) {
						fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totDinucFeatures; cgf++) {
						fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totShapeFeatures; cgf++) {
						fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					//Next, update non-gradient sums
					fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
					rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					//Lastly, add to gradient features
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
									fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
							rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
									rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
						}
					}					
				}
				for (int currOffset=startOffset; currOffset<=endOffset; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					position 		= index-currOffset;
					nucPosOffset	= position*4;
					fDinucPosOffset	= (position-1)*16 + (int) (fNewWorkFeat & 15);
					rDinucPosOffset	= (position-1)*16 + (int) (rNewWorkFeat & 15);
					shapePosOffset 	= (position-2)*nShapeClasses;
					fWeight			= nucAlphas[nucPosOffset + fNewBase];
					rWeight			= nucAlphas[nucPosOffset + rNewBase];
					if (position >= 1) {
						fWeight *= dinucAlphas[fDinucPosOffset];
						rWeight *= dinucAlphas[rDinucPosOffset];
					}
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					//First, update gradient features
					for (int cgf=0; cgf<totNucFeatures; cgf++) {
						fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totDinucFeatures; cgf++) {
						fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totShapeFeatures; cgf++) {
						fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					//Next, update non-gradient sums
					fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
					rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					//Lastly, add to gradient features
					fsNucGradientsWork[nucPosOffset + fNewBase][currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
					rsNucGradientsWork[nucPosOffset + rNewBase][currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					if (position>=1) {
						fsDinucGradientsWork[fDinucPosOffset][currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
						rsDinucGradientsWork[rDinucPosOffset][currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					}
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
									fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
							rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
									rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
						}
					}					
				}
				for (int currOffset=endOffset+1; currOffset<maxFrames; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset 	= (index-currOffset-2)*nShapeClasses;
					fWeight			= 1;
					rWeight			= 1;
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					//First, update gradient features
					for (int cgf=0; cgf<totNucFeatures; cgf++) {
						fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totDinucFeatures; cgf++) {
						fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totShapeFeatures; cgf++) {
						fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					//Next, update non-gradient sums
					fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
					rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					//Lastly, add to gradient features
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
									fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
							rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
									rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
						}
					}					
				}
				
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradients[k][i][j] = fsNucGradientsWork[k][i][j];
						rsNucGradients[k][i][j] = rsNucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradients[k][i][j] = fsDinucGradientsWork[k][i][j];
						rsDinucGradients[k][i][j] = rsDinucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradients[k][i][j] = fsShapeGradientsWork[k][i][j];
						rsShapeGradients[k][i][j] = rsShapeGradientsWork[k][i][j];
					}
				}
			}
		}
		
		//Loop over variable region
		for (int index=flankLength; index<l+flankLength; index++) {								//R0, NUC, DINUC, SHAPE
			startOffset			= (index-k+1 < 0) ? 0 : index-k+1;
			endOffset			= (index > maxFrames-1) ? maxFrames-1 : index;
			shapeStartOffset	= (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset		= (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradientsWork[k][i][j] = 0;
						rsNucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradientsWork[k][i][j] = 0;
						rsDinucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradientsWork[k][i][j] = 0;
						rsShapeGradientsWork[k][i][j] = 0;
					}
				}
			}
			
			for (int currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				for (int newBase=0; newBase<4; newBase++) {
					fNewWorkFeat= (currSubFeat<<2) | newBase;
					rNewWorkFeat= (currSubFeat<<2) | newBase;
					fNextSubFeat= (int) (fNewWorkFeat & subFeatureMask);
					rNextSubFeat= (int) (rNewWorkFeat & subFeatureMask);
					fR0Weight	= r0Alphas[(int) (fNewWorkFeat & r0FeatureMask)];
					rR0Weight	= r0Alphas[(int) reverseComplement((rNewWorkFeat & r0FeatureMask), r0k)];
					for (int currOffset=0; currOffset<startOffset; currOffset++) {
						shapePosOffset 	= (index-currOffset-2)*nShapeClasses;
						fWeight			= fR0Weight;
						rWeight			= rR0Weight;
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							fShapeWeight = 0;
							rShapeWeight = 0;
							fShapeIndex = (int) (fNewWorkFeat & 1023);
							rShapeIndex = (int) (rNewWorkFeat & 1023);
							for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
								fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
								rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
							}
							fWeight *= Math.exp(fShapeWeight);
							rWeight *= Math.exp(rShapeWeight);
						} 
						//First, update gradient features
						for (int cgf=0; cgf<totNucFeatures; cgf++) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						for (int cgf=0; cgf<totDinucFeatures; cgf++) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						for (int cgf=0; cgf<totShapeFeatures; cgf++) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						//Next, update non-gradient sums
						fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
						rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
						//Lastly, add to gradient features
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
								fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
										fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
								rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
										rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
							}
						}
					}
					for (int currOffset=startOffset; currOffset<=endOffset; currOffset++) {
						position		= index-currOffset;
						nucPosOffset	= position*4 + newBase;
						fDinucPosOffset	= (position-1)*16 + (int) (fNewWorkFeat & 15);
						rDinucPosOffset = (position-1)*16 + (int) (rNewWorkFeat & 15);
						shapePosOffset 	= (position-2)*nShapeClasses;
						fWeight			= fR0Weight;
						rWeight			= rR0Weight;
						if (position >= 1) {
							fWeight *= dinucAlphas[fDinucPosOffset];
							rWeight *= dinucAlphas[rDinucPosOffset];								
						}
						fWeight *= nucAlphas[nucPosOffset];
						rWeight *= nucAlphas[nucPosOffset];
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							fShapeWeight = 0;
							rShapeWeight = 0;
							fShapeIndex = (int) (fNewWorkFeat & 1023);
							rShapeIndex = (int) (rNewWorkFeat & 1023);
							for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
								fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
								rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
							}
							fWeight *= Math.exp(fShapeWeight);
							rWeight *= Math.exp(rShapeWeight);
						} 
						//First, update gradient features
						for (int cgf=0; cgf<totNucFeatures; cgf++) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						for (int cgf=0; cgf<totDinucFeatures; cgf++) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						for (int cgf=0; cgf<totShapeFeatures; cgf++) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						//Next, update non-gradient sums
						fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
						rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
						//Lastly, add to gradient features
						fsNucGradientsWork[nucPosOffset][currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;						
						rsNucGradientsWork[nucPosOffset][currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
						if (position>=1) {
							fsDinucGradientsWork[fDinucPosOffset][currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[rDinucPosOffset][currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
						}
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
								fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
										fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
								rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
										rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
							}
						}
					}
					for (int currOffset=endOffset+1; currOffset<maxFrames; currOffset++) {
						shapePosOffset 	= (index-currOffset-2)*nShapeClasses;
						fWeight			= fR0Weight;
						rWeight			= rR0Weight;
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							fShapeWeight = 0;
							rShapeWeight = 0;
							fShapeIndex = (int) (fNewWorkFeat & 1023);
							rShapeIndex = (int) (rNewWorkFeat & 1023);
							for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
								fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
								rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
							}
							fWeight *= Math.exp(fShapeWeight);
							rWeight *= Math.exp(rShapeWeight);
						} 
						//First, update gradient features
						for (int cgf=0; cgf<totNucFeatures; cgf++) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						for (int cgf=0; cgf<totDinucFeatures; cgf++) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						for (int cgf=0; cgf<totShapeFeatures; cgf++) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						//Next, update non-gradient sums
						fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
						rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
						//Lastly, add to gradient features
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
								fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
										fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
								rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
										rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
							}
						}
					}
					
				}
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradients[k][i][j] = fsNucGradientsWork[k][i][j];
						rsNucGradients[k][i][j] = rsNucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradients[k][i][j] = fsDinucGradientsWork[k][i][j];
						rsDinucGradients[k][i][j] = rsDinucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradients[k][i][j] = fsShapeGradientsWork[k][i][j];
						rsShapeGradients[k][i][j] = rsShapeGradientsWork[k][i][j];
					}
				}
			}
		}
		
		//Loop over right fixed region
		for (int index=l+flankLength; index<l+2*flankLength; index++) {								//R0, NUC, DINUC, SHAPE
			startOffset			= (index-k+1 < 0) ? 0 : index-k+1;
			endOffset			= (index > l-k+2*flankLength) ? l-k+2*flankLength : index;
			shapeStartOffset	= (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset		= (index-2 > l-k+2*flankLength) ? l-k+2*flankLength : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradientsWork[k][i][j] = 0;
						rsNucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradientsWork[k][i][j] = 0;
						rsDinucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradientsWork[k][i][j] = 0;
						rsShapeGradientsWork[k][i][j] = 0;
					}
				}
			}
			
			//Find new base
			int fNewBase = (int) fsRFlank & 3;
			int rNewBase = (int) rsRFlank & 3;
			fsRFlank >>= 2;
			rsRFlank >>= 2;
			for (int currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				fNewWorkFeat= (currSubFeat<<2) | fNewBase;
				rNewWorkFeat= (currSubFeat<<2) | rNewBase;
				fNextSubFeat= (int) (fNewWorkFeat & subFeatureMask);
				rNextSubFeat= (int) (rNewWorkFeat & subFeatureMask);
				fShapeIndex	= (int) (fNewWorkFeat & 1023);
				rShapeIndex = (int) (rNewWorkFeat & 1023);
				fR0Weight	= 1;
				rR0Weight	= 1;
				if (index < l+flankLength+r0k-1) {
					fR0Weight = r0Alphas[(int) (fNewWorkFeat & r0FeatureMask)];
					rR0Weight = r0Alphas[(int) reverseComplement((rNewWorkFeat & r0FeatureMask), r0k)];
				}

				for (int currOffset=0; currOffset<startOffset; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset	= (index-currOffset-2)*nShapeClasses;
					fWeight			= fR0Weight;
					rWeight			= rR0Weight;
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					//First, update gradient features
					for (int cgf=0; cgf<totNucFeatures; cgf++) {
						fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totDinucFeatures; cgf++) {
						fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totShapeFeatures; cgf++) {
						fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					//Next, update non-gradient sums
					fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
					rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					//Lastly, add to gradient features
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
									fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
							rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
									rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
						}
					}
				}
				for (int currOffset=startOffset; currOffset<=endOffset; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					position		= index-currOffset;
					nucPosOffset	= position*4;
					fDinucPosOffset = (position-1)*16 + (int) (fNewWorkFeat & 15);
					rDinucPosOffset = (position-1)*16 + (int) (rNewWorkFeat & 15);
					shapePosOffset	= (position-2)*nShapeClasses;
					fWeight			= fR0Weight;
					rWeight			= rR0Weight;
					if (position >= 1) {
						fWeight *= dinucAlphas[fDinucPosOffset];		
						rWeight *= dinucAlphas[rDinucPosOffset];			
					}
					fWeight *= nucAlphas[nucPosOffset + fNewBase];
					rWeight *= nucAlphas[nucPosOffset + rNewBase];
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					//First, update gradient features
					for (int cgf=0; cgf<totNucFeatures; cgf++) {
						fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totDinucFeatures; cgf++) {
						fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totShapeFeatures; cgf++) {
						fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					//Next, update non-gradient sums
					fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
					rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					//Lastly, add to gradient features
					fsNucGradientsWork[nucPosOffset + fNewBase][currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
					rsNucGradientsWork[nucPosOffset + rNewBase][currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					if (position>=1) {
						fsDinucGradientsWork[fDinucPosOffset][currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
						rsDinucGradientsWork[rDinucPosOffset][currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					}
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
									fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
							rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
									rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
						}
					}
				}
				for (int currOffset=endOffset+1; currOffset<maxFrames; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset	= (index-currOffset-2)*nShapeClasses;
					fWeight			= fR0Weight;
					rWeight			= rR0Weight;
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					//First, update gradient features
					for (int cgf=0; cgf<totNucFeatures; cgf++) {
						fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totDinucFeatures; cgf++) {
						fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totShapeFeatures; cgf++) {
						fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					//Next, update non-gradient sums
					fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
					rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					//Lastly, add to gradient features
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
									fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
							rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
									rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
						}
					}
				}
				
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradients[k][i][j] = fsNucGradientsWork[k][i][j];
						rsNucGradients[k][i][j] = rsNucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradients[k][i][j] = fsDinucGradientsWork[k][i][j];
						rsDinucGradients[k][i][j] = rsDinucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradients[k][i][j] = fsShapeGradientsWork[k][i][j];
						rsShapeGradients[k][i][j] = rsShapeGradientsWork[k][i][j];
					}
				}
			}
		}
		
		//Loop over last 2 bases to complete shape features
		int maxIndex = Math.max(2, r0k-1-flankLength);
		for (int index=l+2*flankLength; index<l+2*flankLength+maxIndex; index++) {	
			shapeStartOffset	= (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset		= (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradientsWork[k][i][j] = 0;
						rsNucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradientsWork[k][i][j] = 0;
						rsDinucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradientsWork[k][i][j] = 0;
						rsShapeGradientsWork[k][i][j] = 0;
					}
				}
			}
			
			//Find new base
			int fNewBase = (int) fsRFlank & 3;
			int rNewBase = (int) rsRFlank & 3;
			fsRFlank >>= 2;
			rsRFlank >>= 2;
			for (int currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				fNewWorkFeat= (currSubFeat<<2) | fNewBase;			//What is the new feature?
				rNewWorkFeat= (currSubFeat<<2) | rNewBase;
				fNextSubFeat= (int) (fNewWorkFeat & subFeatureMask);
				rNextSubFeat= (int) (rNewWorkFeat & subFeatureMask);
				fShapeIndex = (int) (fNewWorkFeat & 1023);
				rShapeIndex = (int) (rNewWorkFeat & 1023);
				fR0Weight	= 1;
				rR0Weight	= 1;
				if (index < l+flankLength+r0k-1) {
					fR0Weight = r0Alphas[(int) (fNewWorkFeat & r0FeatureMask)];
					rR0Weight = r0Alphas[(int) reverseComplement((rNewWorkFeat & r0FeatureMask), r0k)];					
				}
				for (int currOffset=0; currOffset<maxFrames; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset	= (index-currOffset-2)*nShapeClasses;
					fWeight			= fR0Weight;
					rWeight			= rR0Weight;
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[fShapeIndex][currShapeFeature];
							rShapeWeight += shapeBetas[shapePosOffset + currShapeFeature]*shapeFeatures[rShapeIndex][currShapeFeature];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					}
					//First, update gradient features
					for (int cgf=0; cgf<totNucFeatures; cgf++) {
						fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totDinucFeatures; cgf++) {
						fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					for (int cgf=0; cgf<totShapeFeatures; cgf++) {
						fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
						rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
					}
					//Next, update non-gradient sums
					fsWork[currOffset][fNextSubFeat] += fsOffsets[currOffset][currSubFeat]*fWeight;
					rsWork[currOffset][rNextSubFeat] += rsOffsets[currOffset][currSubFeat]*rWeight;
					//Lastly, add to gradient features
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int currShapeFeature=0; currShapeFeature<nShapeClasses; currShapeFeature++) {
							fsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][fNextSubFeat] +=
									fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][currShapeFeature];
							rsShapeGradientsWork[shapePosOffset + currShapeFeature][currOffset][rNextSubFeat] +=
									rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][currShapeFeature];
						}
					}
				}
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradients[k][i][j] = fsNucGradientsWork[k][i][j];
						rsNucGradients[k][i][j] = rsNucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradients[k][i][j] = fsDinucGradientsWork[k][i][j];
						rsDinucGradients[k][i][j] = rsDinucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradients[k][i][j] = fsShapeGradientsWork[k][i][j];
						rsShapeGradients[k][i][j] = rsShapeGradientsWork[k][i][j];
					}
				}
			}
		}
		
		Z = 0;
		for (int i=0; i<maxFrames; i++) {
			Z += Array.sum(fsOffsets[i]) + Array.sum(rsOffsets[i]);
		}
		//Sum over frames to get total gradient
		nucGradients = new double[totNucFeatures];
		for (int cgf=0; cgf<totNucFeatures; cgf++) {
			for (int currOffset=0; currOffset<maxFrames; currOffset++) {
				nucGradients[cgf] += Array.sum(fsNucGradients[cgf][currOffset]) 
						+ Array.sum(rsNucGradients[cgf][currOffset]);				
			}
		}
		dinucGradients = new double[totDinucFeatures];
		for (int cgf=0; cgf<totDinucFeatures; cgf++) {
			for (int currOffset=0; currOffset<maxFrames; currOffset++) {
				dinucGradients[cgf] += Array.sum(fsDinucGradients[cgf][currOffset]) 
						+ Array.sum(rsDinucGradients[cgf][currOffset]);				
			}
		}
		shapeGradients = new double[totShapeFeatures];
		for (int cgf=0; cgf<totShapeFeatures; cgf++) {
			for (int currOffset=0; currOffset<maxFrames; currOffset++) {
				shapeGradients[cgf] += Array.sum(fsShapeGradients[cgf][currOffset])
						+ Array.sum(rsShapeGradients[cgf][currOffset]);
			}
		}
		fsOffsets			= null;
		rsOffsets			= null;
		fsWork				= null;
		rsWork				= null;
		fsNucGradients		= null;
		rsNucGradients		= null;
		fsDinucGradients	= null;
		rsDinucGradients	= null;
		fsShapeGradients 	= null;
		rsShapeGradients	= null;
		fsNucGradientsWork	= null;
		rsNucGradientsWork	= null;
		fsDinucGradientsWork= null;
		rsDinucGradientsWork= null;
		fsShapeGradientsWork= null;
		rsShapeGradientsWork= null;
		System.gc();
	}

	@Override
	public void recursiveHessian() {
		int startOffset;
		int endOffset;
		int shapeStartOffset;
		int shapeEndOffset;
		int position;
		int nucPosOffset, fNucPosOffset, rNucPosOffset;
		int fDinucPosOffset;
		int rDinucPosOffset;
		int shapePosOffset;
		int fNextSubFeat;
		int rNextSubFeat;
		int fShapeIndex					= 0;
		int rShapeIndex					= 0;
		int tempIdx;
		double fwdUpdate, revUpdate;
		double fWeight					= 0;
		double rWeight					= 0;
		double fShapeWeight				= 0;
		double rShapeWeight				= 0;
		double fR0Weight				= 0;
		double rR0Weight				= 0;
		double[][] fsOffsets			= matrixClone(fsInitMatrix);
		double[][] rsOffsets			= matrixClone(rsInitMatrix);
		double[][] fsWork				= new double[maxFrames][maxSubFeatures];
		double[][] rsWork				= new double[maxFrames][maxSubFeatures];
		double[][][] fsNucGradients		= new double[totNucFeatures][maxFrames][maxSubFeatures];
		double[][][] rsNucGradients		= new double[totNucFeatures][maxFrames][maxSubFeatures];
		double[][][] fsDinucGradients	= new double[totDinucFeatures][maxFrames][maxSubFeatures];
		double[][][] rsDinucGradients	= new double[totDinucFeatures][maxFrames][maxSubFeatures];
		double[][][] fsShapeGradients	= new double[totShapeFeatures][maxFrames][maxSubFeatures];
		double[][][] rsShapeGradients	= new double[totShapeFeatures][maxFrames][maxSubFeatures];
		double[][][] fsNucGradientsWork = new double[totNucFeatures][maxFrames][maxSubFeatures];
		double[][][] rsNucGradientsWork = new double[totNucFeatures][maxFrames][maxSubFeatures];
		double[][][] fsDinucGradientsWork=new double[totDinucFeatures][maxFrames][maxSubFeatures];
		double[][][] rsDinucGradientsWork=new double[totDinucFeatures][maxFrames][maxSubFeatures];
		double[][][] fsShapeGradientsWork=new double[totShapeFeatures][maxFrames][maxSubFeatures];
		double[][][] rsShapeGradientsWork=new double[totShapeFeatures][maxFrames][maxSubFeatures];
		double[][][][] fsHessian		= new double[totFeatures][totFeatures][maxFrames][maxSubFeatures];
		double[][][][] rsHessian		= new double[totFeatures][totFeatures][maxFrames][maxSubFeatures];
		double[][][][] fsHessianWork	= new double[totFeatures][totFeatures][maxFrames][maxSubFeatures];
		double[][][][] rsHessianWork	= new double[totFeatures][totFeatures][maxFrames][maxSubFeatures];
		long fNewWorkFeat;
		long rNewWorkFeat;
		long fsLFlank					= leftFlank;
		long fsRFlank					= rightFlank;
		long rsLFlank					= leftFlankRC;
		long rsRFlank					= rightFlankRC;
		
		//Loop over left fixed region
		for (int index=0; index<flankLength; index++) {									//NO ROUND0 MODEL HERE
			startOffset			= (index-k+1 < 0) ? 0 : index-k+1;
			endOffset			= (index > maxFrames-1) ? maxFrames-1 : index;
			shapeStartOffset	= (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset		= (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradientsWork[k][i][j] = 0;
						rsNucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradientsWork[k][i][j] = 0;
						rsDinucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradientsWork[k][i][j] = 0;
						rsShapeGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int l=0; l<totFeatures; l++) {
				for (int k=0; k<totFeatures; k++) {
					for (int i=0; i<maxFrames; i++) {
						for (int j=0; j<maxSubFeatures; j++) {
							fsHessianWork[l][k][i][j] = 0;
							rsHessianWork[l][k][i][j] = 0;
						}
					}
				}				
			}
			
			//Find new base
			int fNewBase = (int) fsLFlank & 3;
			int rNewBase = (int) rsLFlank & 3;
			fsLFlank >>= 2;
			rsLFlank >>= 2;
			for (int currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				fNewWorkFeat= (currSubFeat<<2) | fNewBase;
				rNewWorkFeat= (currSubFeat<<2) | rNewBase;
				fNextSubFeat= (int) (fNewWorkFeat & subFeatureMask);
				rNextSubFeat= (int) (rNewWorkFeat & subFeatureMask);
				fShapeIndex	= (int) (fNewWorkFeat & 1023);
				rShapeIndex = (int) (rNewWorkFeat & 1023);
				
				for (int currOffset=0; currOffset<startOffset; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset 	= (index-currOffset-2)*nShapeClasses;
					fWeight			= 1;
					rWeight			= 1;
					//Adjust weights if shape features exist
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int csc=0; csc<nShapeClasses; csc++) {
							fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
							rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					//Update all positions
					for (int cgf=0; cgf<totFeatures; cgf++) {
						if (cgf<totNucFeatures) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totDinucFeatures) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totShapeFeatures) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
						}
						for (int cgf2=0; cgf2<=cgf; cgf2++) {
							fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
							rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
						}
					}
					fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
					revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
					fsWork[currOffset][fNextSubFeat] += fwdUpdate;
					rsWork[currOffset][rNextSubFeat] += revUpdate;
					
					//If shape features exist, add appropriate gradient components
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int csc=0; csc<nShapeClasses; csc++) {
							fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
									fwdUpdate*shapeFeatures[fShapeIndex][csc];
							rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
									revUpdate*shapeFeatures[rShapeIndex][csc];	
						}
						shapePosOffset += 4*k + 16*(k-1);
						for (int csc1=0; csc1<nShapeClasses; csc1++) {
							for (int csc2=0; csc2<=csc1; csc2++) {
								fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
								rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
							}
						}
						
						for (int csc=0; csc<nShapeClasses; csc++) {
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[shapePosOffset + csc][cgf][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][cgf][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf + 4*k;
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
						}
					}					
				}
				//Active window.
				for (int currOffset=startOffset; currOffset<=endOffset; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					position 		= index-currOffset;
					fNucPosOffset	= position*4 + fNewBase;
					rNucPosOffset	= position*4 + rNewBase;
					fDinucPosOffset	= (position-1)*16 + (int) (fNewWorkFeat & 15);
					rDinucPosOffset	= (position-1)*16 + (int) (rNewWorkFeat & 15);
					shapePosOffset 	= (position-2)*nShapeClasses;
					fWeight			= nucAlphas[fNucPosOffset];
					rWeight			= nucAlphas[rNucPosOffset];
					if (position >= 1) {
						fWeight *= dinucAlphas[fDinucPosOffset];
						rWeight *= dinucAlphas[rDinucPosOffset];
					}
					//Adjust weights if shape features exist
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int csc=0; csc<nShapeClasses; csc++) {
							fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
							rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					}
					
					//Update all positions
					for (int cgf=0; cgf<totFeatures; cgf++) {
						if (cgf<totNucFeatures) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totDinucFeatures) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totShapeFeatures) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
						}
						for (int cgf2=0; cgf2<=cgf; cgf2++) {
							fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
							rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
						}
					}
					fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
					revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
					fsWork[currOffset][fNextSubFeat]									+= fwdUpdate;
					rsWork[currOffset][rNextSubFeat]									+= revUpdate;						
					//Add to Nucleotide Gradients
					fsNucGradientsWork[fNucPosOffset][currOffset][fNextSubFeat]			+= fwdUpdate;						
					rsNucGradientsWork[rNucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
					//Add to Nucleotide Hessians
					fsHessianWork[fNucPosOffset][fNucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
					rsHessianWork[rNucPosOffset][rNucPosOffset][currOffset][rNextSubFeat] += revUpdate;
										
					//Do Shape features exist?
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						//Add to Shape Gradients
						for (int csc=0; csc<nShapeClasses; csc++) {
							fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
									fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][csc];
							rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
									rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][csc];
						}
						//Add to Shape Hessians
						shapePosOffset += 4*k + 16*(k-1);
						for (int csc1=0; csc1<nShapeClasses; csc1++) {
							for (int csc2=0; csc2<=csc1; csc2++) {
								fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
								rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
							}
							//Nuc/Shape Hessian Cross Terms
							fsHessianWork[shapePosOffset + csc1][fNucPosOffset][currOffset][fNextSubFeat] +=
									fwdUpdate*shapeFeatures[fShapeIndex][csc1];
							rsHessianWork[shapePosOffset + csc1][rNucPosOffset][currOffset][rNextSubFeat] +=
									revUpdate*shapeFeatures[rShapeIndex][csc1];
						}
						
						//Do dinucleotides exist?
						if (position>=1) {
							//Add to Dinucleotide Gradients
							fsDinucGradientsWork[fDinucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
							rsDinucGradientsWork[rDinucPosOffset][currOffset][rNextSubFeat] += revUpdate;
							//Add to Dinucleotide Hessians
							fDinucPosOffset += 4*k;
							rDinucPosOffset += 4*k;
							fsHessianWork[fDinucPosOffset][fDinucPosOffset][currOffset][fNextSubFeat]	+= fwdUpdate;
							rsHessianWork[rDinucPosOffset][rDinucPosOffset][currOffset][rNextSubFeat]	+= revUpdate;
							//Add to Nucleotide/Dinucleotide Hessian Cross Terms
							fsHessianWork[fDinucPosOffset][fNucPosOffset][currOffset][fNextSubFeat]		+= fwdUpdate;
							rsHessianWork[rDinucPosOffset][rNucPosOffset][currOffset][rNextSubFeat]		+= revUpdate;
							//Add to Shape/Dinucleotide Hessian Cross Terms
							for (int csc=0; csc<nShapeClasses; csc++) {
								fsHessianWork[shapePosOffset + csc][fDinucPosOffset][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc];
								rsHessianWork[shapePosOffset + csc][rDinucPosOffset][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc];
							}
							//Carry over all values from previous gradient
							//Carry over nucleotide features from all previous features into the nucleotide bin with the current features
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[fNucPosOffset][cgf][currOffset][fNextSubFeat]				+= fwdUpdate;
								rsHessianWork[rNucPosOffset][cgf][currOffset][rNextSubFeat]				+= revUpdate;
								fsHessianWork[fDinucPosOffset][cgf][currOffset][fNextSubFeat]			+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][cgf][currOffset][rNextSubFeat]			+= revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][cgf][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][cgf][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}
							}
							//Same, but with dinucleotide features
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k;
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat]			+= fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
								fsHessianWork[fDinucPosOffset][tempIdx][currOffset][fNextSubFeat]		+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][tempIdx][currOffset][rNextSubFeat]		+= revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}
							}
							//Same, but with shape features
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat] 			+= fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
								fsHessianWork[tempIdx][fDinucPosOffset][currOffset][fNextSubFeat]		+= fwdUpdate;
								rsHessianWork[tempIdx][rDinucPosOffset][currOffset][rNextSubFeat]		+= revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}		
							}
						} else {	//No Dinucleotides exist.
							//Carry over all values from previous gradient
							//Carry over nucleotide features from all previous features into the nucleotide bin with the current features 
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[rNucPosOffset][cgf][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[fNucPosOffset][cgf][currOffset][rNextSubFeat] += revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][cgf][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][cgf][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}
							}
							//Same, but with dinucleotide features
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k;
								fsHessianWork[fNucPosOffset][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[rNucPosOffset][tempIdx][currOffset][rNextSubFeat] += revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}
							}
							//Same, but with shape features
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat] += revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}		
							}
						}
					} else {	//No Shape features exist.
						//Do dinucleotides exist?
						if (position>=1) {
							//Add to Dinucleotide Gradients
							fsDinucGradientsWork[fDinucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
							rsDinucGradientsWork[rDinucPosOffset][currOffset][rNextSubFeat] += revUpdate;
							//Add to Dinucleotide Hessians
							fDinucPosOffset += 4*k;
							rDinucPosOffset += 4*k;
							fsHessianWork[fDinucPosOffset][fDinucPosOffset][currOffset][fNextSubFeat]	+= fwdUpdate;
							rsHessianWork[rDinucPosOffset][rDinucPosOffset][currOffset][rNextSubFeat]	+= revUpdate;
							//Add to Nucleotide/Dinucleotide Hessian Cross Terms
							fsHessianWork[fDinucPosOffset][fNucPosOffset][currOffset][fNextSubFeat]		+= fwdUpdate;
							rsHessianWork[rDinucPosOffset][rNucPosOffset][currOffset][rNextSubFeat]		+= revUpdate;
							//Carry over all values from previous gradient
							//Carry over nucleotide features from all previous features into the nucleotide bin with the current features
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[fNucPosOffset][cgf][currOffset][fNextSubFeat]				+= fwdUpdate;
								rsHessianWork[rNucPosOffset][cgf][currOffset][rNextSubFeat]				+= revUpdate;
								fsHessianWork[fDinucPosOffset][cgf][currOffset][fNextSubFeat]			+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][cgf][currOffset][rNextSubFeat]			+= revUpdate;
							}
							//Same, but with dinucleotide features
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k;
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat]			+= fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
								fsHessianWork[fDinucPosOffset][tempIdx][currOffset][fNextSubFeat]		+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][tempIdx][currOffset][rNextSubFeat]		+= revUpdate;
							}
							//Same, but with shape features
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat] += revUpdate;
								fsHessianWork[tempIdx][fDinucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[tempIdx][rDinucPosOffset][currOffset][rNextSubFeat] += revUpdate;
							}
						} else {	//No Dinucleotides exist.
							//Carry over all values from previous gradient
							//Carry over nucleotide features from all previous features into the nucleotide bin with the current features 
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[fNucPosOffset][cgf][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[rNucPosOffset][cgf][currOffset][rNextSubFeat] += revUpdate;
							}
							//Same, but with dinucleotide features
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k;
								fsHessianWork[fNucPosOffset][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[rNucPosOffset][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat] += revUpdate;
							}
						}
					}					
				}
				
				//Not in an active window.
				for (int currOffset=endOffset+1; currOffset<maxFrames; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset 	= (index-currOffset-2)*nShapeClasses;
					fWeight			= 1;
					rWeight			= 1;
					//Adjust weights if shape features exist
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int csc=0; csc<nShapeClasses; csc++) {
							fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
							rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					}
					//Update all positions
					for (int cgf=0; cgf<totFeatures; cgf++) {
						if (cgf<totNucFeatures) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totDinucFeatures) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totShapeFeatures) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
						}
						for (int cgf2=0; cgf2<=cgf; cgf2++) {
							fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
							rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
						}
					}
					fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
					revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
					fsWork[currOffset][fNextSubFeat] += fwdUpdate;
					rsWork[currOffset][rNextSubFeat] += revUpdate;
					
					//If shape features exist, add appropriate gradient components
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int csc=0; csc<nShapeClasses; csc++) {
							fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
									fwdUpdate*shapeFeatures[fShapeIndex][csc];
							rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
									revUpdate*shapeFeatures[rShapeIndex][csc];	
						}
						shapePosOffset += 4*k + 16*(k-1);
						for (int csc1=0; csc1<nShapeClasses; csc1++) {
							for (int csc2=0; csc2<=csc1; csc2++) {
								fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
								rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
							}
						}
						
						for (int csc=0; csc<nShapeClasses; csc++) {
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[shapePosOffset + csc][cgf][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][cgf][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf + 4*k;
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
						}
					}					
				}
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradients[k][i][j] = fsNucGradientsWork[k][i][j];
						rsNucGradients[k][i][j] = rsNucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradients[k][i][j] = fsDinucGradientsWork[k][i][j];
						rsDinucGradients[k][i][j] = rsDinucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradients[k][i][j] = fsShapeGradientsWork[k][i][j];
						rsShapeGradients[k][i][j] = rsShapeGradientsWork[k][i][j];
					}
				}
			}
			for (int l=0; l<totFeatures; l++) {
				for (int k=0; k<totFeatures; k++) {
					for (int i=0; i<maxFrames; i++) {
						for (int j=0; j<maxSubFeatures; j++) {
							fsHessian[l][k][i][j] = fsHessianWork[l][k][i][j];
							rsHessian[l][k][i][j] = rsHessianWork[l][k][i][j];
						}
					}
				}				
			}
		}
		
		//Loop over variable region
		for (int index=flankLength; index<l+flankLength; index++) {								//R0, NUC, DINUC, SHAPE
			startOffset			= (index-k+1 < 0) ? 0 : index-k+1;
			endOffset			= (index > maxFrames-1) ? maxFrames-1 : index;
			shapeStartOffset	= (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset		= (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradientsWork[k][i][j] = 0;
						rsNucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradientsWork[k][i][j] = 0;
						rsDinucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradientsWork[k][i][j] = 0;
						rsShapeGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int l=0; l<totFeatures; l++) {
				for (int k=0; k<totFeatures; k++) {
					for (int i=0; i<maxFrames; i++) {
						for (int j=0; j<maxSubFeatures; j++) {
							fsHessianWork[l][k][i][j] = 0;
							rsHessianWork[l][k][i][j] = 0;
						}
					}
				}				
			}
			
			for (int currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				for (int newBase=0; newBase<4; newBase++) {
					fNewWorkFeat= (currSubFeat<<2) | newBase;
					rNewWorkFeat= (currSubFeat<<2) | newBase;
					fNextSubFeat= (int) (fNewWorkFeat & subFeatureMask);
					rNextSubFeat= (int) (rNewWorkFeat & subFeatureMask);
					fShapeIndex = (int) (fNewWorkFeat & 1023);
					rShapeIndex = (int) (rNewWorkFeat & 1023);
					fR0Weight	= r0Alphas[(int) (fNewWorkFeat & r0FeatureMask)];
					rR0Weight	= r0Alphas[(int) reverseComplement((rNewWorkFeat & r0FeatureMask), r0k)];
					//Not in an active window.
					for (int currOffset=0; currOffset<startOffset; currOffset++) {
						if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
						shapePosOffset 	= (index-currOffset-2)*nShapeClasses;
						fWeight			= fR0Weight;
						rWeight			= rR0Weight;
						//Adjust weights if shape features exist
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							fShapeWeight = 0;
							rShapeWeight = 0;
							for (int csc=0; csc<nShapeClasses; csc++) {
								fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
								rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
							}
							fWeight *= Math.exp(fShapeWeight);
							rWeight *= Math.exp(rShapeWeight);
						}
						//Update all positions
						for (int cgf=0; cgf<totFeatures; cgf++) {
							if (cgf<totNucFeatures) {
								fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
							}
							if (cgf<totDinucFeatures) {
								fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
							}
							if (cgf<totShapeFeatures) {
								fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
							}
							for (int cgf2=0; cgf2<=cgf; cgf2++) {
								fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
								rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
							}
						}
						fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
						revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
						fsWork[currOffset][fNextSubFeat] += fwdUpdate;
						rsWork[currOffset][rNextSubFeat] += revUpdate;
						
						//If shape features exist, add appropriate gradient components
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							for (int csc=0; csc<nShapeClasses; csc++) {
								fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc];
								rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc];	
							}
							shapePosOffset += 4*k + 16*(k-1);
							for (int csc1=0; csc1<nShapeClasses; csc1++) {
								for (int csc2=0; csc2<=csc1; csc2++) {
									fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
									rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
								}
							}
							
							for (int csc=0; csc<nShapeClasses; csc++) {
								for (int cgf=0; cgf<totNucFeatures; cgf++) {
									fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = shapeFeatures[rShapeIndex][csc]*rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									fsHessianWork[shapePosOffset + csc][cgf][currOffset][fNextSubFeat]	+= fwdUpdate;
									rsHessianWork[shapePosOffset + csc][cgf][currOffset][rNextSubFeat] += revUpdate;
								}
								for (int cgf=0; cgf<totDinucFeatures; cgf++) {
									fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = shapeFeatures[rShapeIndex][csc]*rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf + 4*k;
									fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat]	+= fwdUpdate;
									rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
								}
								for (int cgf=0; cgf<totShapeFeatures; cgf++) {
									fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = shapeFeatures[rShapeIndex][csc]*rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k+16*(k-1);
									fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
									rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
								}
							}
						}
					}
					//Active window.
					for (int currOffset=startOffset; currOffset<=endOffset; currOffset++) {
						position		= index-currOffset;
						nucPosOffset	= position*4 + newBase;
						fDinucPosOffset	= (position-1)*16 + (int) (fNewWorkFeat & 15);
						rDinucPosOffset = (position-1)*16 + (int) (rNewWorkFeat & 15);
						shapePosOffset 	= (position-2)*nShapeClasses;
						fWeight			= fR0Weight;
						rWeight			= rR0Weight;
						if (position >= 1) {
							fWeight *= dinucAlphas[fDinucPosOffset];
							rWeight *= dinucAlphas[rDinucPosOffset];								
						}
						fWeight *= nucAlphas[nucPosOffset];
						rWeight *= nucAlphas[nucPosOffset];
						//Adjust weights if shape features exist
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							fShapeWeight = 0;
							rShapeWeight = 0;
							for (int csc=0; csc<nShapeClasses; csc++) {
								fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
								rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
							}
							fWeight *= Math.exp(fShapeWeight);
							rWeight *= Math.exp(rShapeWeight);
						}
						
						//Update all positions
						for (int cgf=0; cgf<totFeatures; cgf++) {
							if (cgf<totNucFeatures) {
								fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
							}
							if (cgf<totDinucFeatures) {
								fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
							}
							if (cgf<totShapeFeatures) {
								fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
							}
							for (int cgf2=0; cgf2<=cgf; cgf2++) {
								fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
								rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
							}
						}
						fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
						revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
						fsWork[currOffset][fNextSubFeat]									+= fwdUpdate;
						rsWork[currOffset][rNextSubFeat]									+= revUpdate;						
						//Add to Nucleotide Gradients
						fsNucGradientsWork[nucPosOffset][currOffset][fNextSubFeat]			+= fwdUpdate;						
						rsNucGradientsWork[nucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
						//Add to Nucleotide Hessians
						fsHessianWork[nucPosOffset][nucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
						rsHessianWork[nucPosOffset][nucPosOffset][currOffset][rNextSubFeat] += revUpdate;
						
						//Do Shape features exist?
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							//Add to Shape Gradients
							for (int csc=0; csc<nShapeClasses; csc++) {
								fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
										fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][csc];
								rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
										rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][csc];
							}
							//Add to Shape Hessians
							shapePosOffset += 4*k + 16*(k-1);
							for (int csc1=0; csc1<nShapeClasses; csc1++) {
								for (int csc2=0; csc2<=csc1; csc2++) {
									fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
									rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
								}
								//Nuc/Shape Hessian Cross Terms
								fsHessianWork[shapePosOffset + csc1][nucPosOffset][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc1];
								rsHessianWork[shapePosOffset + csc1][nucPosOffset][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc1];
							}
							
							//Do dinucleotides exist?
							if (position>=1) {
								//Add to Dinucleotide Gradients
								fsDinucGradientsWork[fDinucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
								rsDinucGradientsWork[rDinucPosOffset][currOffset][rNextSubFeat] += revUpdate;
								//Add to Dinucleotide Hessians
								fDinucPosOffset += 4*k;
								rDinucPosOffset += 4*k;
								fsHessianWork[fDinucPosOffset][fDinucPosOffset][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][rDinucPosOffset][currOffset][rNextSubFeat]	+= revUpdate;
								//Add to Nucleotide/Dinucleotide Hessian Cross Terms
								fsHessianWork[fDinucPosOffset][nucPosOffset][currOffset][fNextSubFeat]		+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][nucPosOffset][currOffset][rNextSubFeat]		+= revUpdate;
								//Add to Shape/Dinucleotide Hessian Cross Terms
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset + csc][fDinucPosOffset][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset + csc][rDinucPosOffset][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}
								//Carry over all values from previous gradient
								//Carry over nucleotide features from all previous features into the nucleotide bin with the current features
								for (int cgf=0; cgf<totNucFeatures; cgf++) {
									fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									fsHessianWork[nucPosOffset][cgf][currOffset][fNextSubFeat]				+= fwdUpdate;
									rsHessianWork[nucPosOffset][cgf][currOffset][rNextSubFeat]				+= revUpdate;
									fsHessianWork[fDinucPosOffset][cgf][currOffset][fNextSubFeat]			+= fwdUpdate;
									rsHessianWork[rDinucPosOffset][cgf][currOffset][rNextSubFeat]			+= revUpdate;
									for (int csc=0; csc<nShapeClasses; csc++) {
										fsHessianWork[shapePosOffset+csc][cgf][currOffset][fNextSubFeat] +=
												fwdUpdate*shapeFeatures[fShapeIndex][csc];
										rsHessianWork[shapePosOffset+csc][cgf][currOffset][rNextSubFeat] +=
												revUpdate*shapeFeatures[rShapeIndex][csc];
									}
								}
								//Same, but with dinucleotide features
								for (int cgf=0; cgf<totDinucFeatures; cgf++) {
									fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k;
									fsHessianWork[tempIdx][nucPosOffset][currOffset][fNextSubFeat]			+= fwdUpdate;
									rsHessianWork[tempIdx][nucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
									fsHessianWork[fDinucPosOffset][tempIdx][currOffset][fNextSubFeat]		+= fwdUpdate;
									rsHessianWork[rDinucPosOffset][tempIdx][currOffset][rNextSubFeat]		+= revUpdate;
									for (int csc=0; csc<nShapeClasses; csc++) {
										fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
												fwdUpdate*shapeFeatures[fShapeIndex][csc];
										rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
												revUpdate*shapeFeatures[rShapeIndex][csc];
									}
								}
								//Same, but with shape features
								for (int cgf=0; cgf<totShapeFeatures; cgf++) {
									fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k+16*(k-1);
									fsHessianWork[tempIdx][nucPosOffset][currOffset][fNextSubFeat] 			+= fwdUpdate;
									rsHessianWork[tempIdx][nucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
									fsHessianWork[tempIdx][fDinucPosOffset][currOffset][fNextSubFeat]		+= fwdUpdate;
									rsHessianWork[tempIdx][rDinucPosOffset][currOffset][rNextSubFeat]		+= revUpdate;
									for (int csc=0; csc<nShapeClasses; csc++) {
										fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
												fwdUpdate*shapeFeatures[fShapeIndex][csc];
										rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
												revUpdate*shapeFeatures[rShapeIndex][csc];
									}		
								}
							} else {	//No Dinucleotides exist.
								//Carry over all values from previous gradient
								//Carry over nucleotide features from all previous features into the nucleotide bin with the current features 
								for (int cgf=0; cgf<totNucFeatures; cgf++) {
									fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									fsHessianWork[nucPosOffset][cgf][currOffset][fNextSubFeat] += fwdUpdate;
									rsHessianWork[nucPosOffset][cgf][currOffset][rNextSubFeat] += revUpdate;
									for (int csc=0; csc<nShapeClasses; csc++) {
										fsHessianWork[shapePosOffset+csc][cgf][currOffset][fNextSubFeat] +=
												fwdUpdate*shapeFeatures[fShapeIndex][csc];
										rsHessianWork[shapePosOffset+csc][cgf][currOffset][rNextSubFeat] +=
												revUpdate*shapeFeatures[rShapeIndex][csc];
									}
								}
								//Same, but with dinucleotide features
								for (int cgf=0; cgf<totDinucFeatures; cgf++) {
									fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k;
									fsHessianWork[nucPosOffset][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
									rsHessianWork[nucPosOffset][tempIdx][currOffset][rNextSubFeat] += revUpdate;
									for (int csc=0; csc<nShapeClasses; csc++) {
										fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
												fwdUpdate*shapeFeatures[fShapeIndex][csc];
										rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
												revUpdate*shapeFeatures[rShapeIndex][csc];
									}
								}
								//Same, but with shape features
								for (int cgf=0; cgf<totShapeFeatures; cgf++) {
									fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k+16*(k-1);
									fsHessianWork[tempIdx][nucPosOffset][currOffset][fNextSubFeat] +=
											fwdUpdate;//*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[tempIdx][nucPosOffset][currOffset][rNextSubFeat] +=
											revUpdate;//*shapeFeatures[rShapeIndex][csc];
									for (int csc=0; csc<nShapeClasses; csc++) {
										fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
												fwdUpdate*shapeFeatures[fShapeIndex][csc];
										rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
												revUpdate*shapeFeatures[rShapeIndex][csc];
									}		
								}
							}
						} else {	//No Shape features exist.
							//Do dinucleotides exist?
							if (position>=1) {
								//Add to Dinucleotide Gradients
								fsDinucGradientsWork[fDinucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
								rsDinucGradientsWork[rDinucPosOffset][currOffset][rNextSubFeat] += revUpdate;
								//Add to Dinucleotide Hessians
								fDinucPosOffset += 4*k;
								rDinucPosOffset += 4*k;
								fsHessianWork[fDinucPosOffset][fDinucPosOffset][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][rDinucPosOffset][currOffset][rNextSubFeat]	+= revUpdate;
								//Add to Nucleotide/Dinucleotide Hessian Cross Terms
								fsHessianWork[fDinucPosOffset][nucPosOffset][currOffset][fNextSubFeat]		+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][nucPosOffset][currOffset][rNextSubFeat]		+= revUpdate;
								//Carry over all values from previous gradient
								//Carry over nucleotide features from all previous features into the nucleotide bin with the current features
								for (int cgf=0; cgf<totNucFeatures; cgf++) {
									fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									fsHessianWork[nucPosOffset][cgf][currOffset][fNextSubFeat]				+= fwdUpdate;
									rsHessianWork[nucPosOffset][cgf][currOffset][rNextSubFeat]				+= revUpdate;
									fsHessianWork[fDinucPosOffset][cgf][currOffset][fNextSubFeat]			+= fwdUpdate;
									rsHessianWork[rDinucPosOffset][cgf][currOffset][rNextSubFeat]			+= revUpdate;
								}
								//Same, but with dinucleotide features
								for (int cgf=0; cgf<totDinucFeatures; cgf++) {
									fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k;
									fsHessianWork[tempIdx][nucPosOffset][currOffset][fNextSubFeat]			+= fwdUpdate;
									rsHessianWork[tempIdx][nucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
									fsHessianWork[fDinucPosOffset][tempIdx][currOffset][fNextSubFeat]		+= fwdUpdate;
									rsHessianWork[rDinucPosOffset][tempIdx][currOffset][rNextSubFeat]		+= revUpdate;
								}
								//Same, but with shape features
								for (int cgf=0; cgf<totShapeFeatures; cgf++) {
									fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k+16*(k-1);
									fsHessianWork[tempIdx][nucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
									rsHessianWork[tempIdx][nucPosOffset][currOffset][rNextSubFeat] += revUpdate;
									fsHessianWork[tempIdx][fDinucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
									rsHessianWork[tempIdx][rDinucPosOffset][currOffset][rNextSubFeat] += revUpdate;
								}
							} else {	//No Dinucleotides exist.
								//Carry over all values from previous gradient
								//Carry over nucleotide features from all previous features into the nucleotide bin with the current features 
								for (int cgf=0; cgf<totNucFeatures; cgf++) {
									fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									fsHessianWork[nucPosOffset][cgf][currOffset][fNextSubFeat] += fwdUpdate;
									rsHessianWork[nucPosOffset][cgf][currOffset][rNextSubFeat] += revUpdate;
								}
								//Same, but with dinucleotide features
								for (int cgf=0; cgf<totDinucFeatures; cgf++) {
									fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k;
									fsHessianWork[nucPosOffset][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
									rsHessianWork[nucPosOffset][tempIdx][currOffset][rNextSubFeat] += revUpdate;
								}
								for (int cgf=0; cgf<totShapeFeatures; cgf++) {
									fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k+16*(k-1);
									fsHessianWork[tempIdx][nucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
									rsHessianWork[tempIdx][nucPosOffset][currOffset][rNextSubFeat] += revUpdate;
								}
							}
						}
					}

					//Not in an active window.
					for (int currOffset=endOffset+1; currOffset<maxFrames; currOffset++) {
						shapePosOffset 	= (index-currOffset-2)*nShapeClasses;
						fWeight			= fR0Weight;
						rWeight			= rR0Weight;
						//Adjust weights if shape features exist
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							fShapeWeight = 0;
							rShapeWeight = 0;
							for (int csc=0; csc<nShapeClasses; csc++) {
								fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
								rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
							}
							fWeight *= Math.exp(fShapeWeight);
							rWeight *= Math.exp(rShapeWeight);
						}
						//Update all positions
						for (int cgf=0; cgf<totFeatures; cgf++) {
							if (cgf<totNucFeatures) {
								fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
							}
							if (cgf<totDinucFeatures) {
								fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
							}
							if (cgf<totShapeFeatures) {
								fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
							}
							for (int cgf2=0; cgf2<=cgf; cgf2++) {
								fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
								rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
							}
						}
						fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
						revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
						fsWork[currOffset][fNextSubFeat] += fwdUpdate;
						rsWork[currOffset][rNextSubFeat] += revUpdate;
						//If shape features exist, add appropriate gradient components
						if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
							for (int csc=0; csc<nShapeClasses; csc++) {
								fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc];
								rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc];	
							}
							shapePosOffset += 4*k + 16*(k-1);
							for (int csc1=0; csc1<nShapeClasses; csc1++) {
								for (int csc2=0; csc2<=csc1; csc2++) {
									fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
									rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
								}
							}
							
							for (int csc=0; csc<nShapeClasses; csc++) {
								for (int cgf=0; cgf<totNucFeatures; cgf++) {
									fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = shapeFeatures[rShapeIndex][csc]*rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									fsHessianWork[shapePosOffset + csc][cgf][currOffset][fNextSubFeat]	+= fwdUpdate;
									rsHessianWork[shapePosOffset + csc][cgf][currOffset][rNextSubFeat] += revUpdate;
								}
								for (int cgf=0; cgf<totDinucFeatures; cgf++) {
									fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = shapeFeatures[rShapeIndex][csc]*rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf + 4*k;
									fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat]	+= fwdUpdate;
									rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
								}
								for (int cgf=0; cgf<totShapeFeatures; cgf++) {
									fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
									revUpdate = shapeFeatures[rShapeIndex][csc]*rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
									if (fwdUpdate==0 && revUpdate==0)	continue;
									tempIdx = cgf+4*k+16*(k-1);
									fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
									rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
								}
							}
						}
					}
				}
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradients[k][i][j] = fsNucGradientsWork[k][i][j];
						rsNucGradients[k][i][j] = rsNucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradients[k][i][j] = fsDinucGradientsWork[k][i][j];
						rsDinucGradients[k][i][j] = rsDinucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradients[k][i][j] = fsShapeGradientsWork[k][i][j];
						rsShapeGradients[k][i][j] = rsShapeGradientsWork[k][i][j];
					}
				}
			}
			for (int l=0; l<totFeatures; l++) {
				for (int k=0; k<totFeatures; k++) {
					for (int i=0; i<maxFrames; i++) {
						for (int j=0; j<maxSubFeatures; j++) {
							fsHessian[l][k][i][j] = fsHessianWork[l][k][i][j];
							rsHessian[l][k][i][j] = rsHessianWork[l][k][i][j];
						}
					}
				}				
			}
		}
		
		//Loop over right fixed region
		for (int index=l+flankLength; index<l+2*flankLength; index++) {								//R0, NUC, DINUC, SHAPE
			startOffset			= (index-k+1 < 0) ? 0 : index-k+1;
			endOffset			= (index > l-k+2*flankLength) ? l-k+2*flankLength : index;
			shapeStartOffset	= (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset		= (index-2 > l-k+2*flankLength) ? l-k+2*flankLength : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradientsWork[k][i][j] = 0;
						rsNucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradientsWork[k][i][j] = 0;
						rsDinucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradientsWork[k][i][j] = 0;
						rsShapeGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int l=0; l<totFeatures; l++) {
				for (int k=0; k<totFeatures; k++) {
					for (int i=0; i<maxFrames; i++) {
						for (int j=0; j<maxSubFeatures; j++) {
							fsHessianWork[l][k][i][j] = 0;
							rsHessianWork[l][k][i][j] = 0;
						}
					}
				}				
			}
			
			//Find new base
			int fNewBase = (int) fsRFlank & 3;
			int rNewBase = (int) rsRFlank & 3;
			fsRFlank >>= 2;
			rsRFlank >>= 2;
			for (int currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				fNewWorkFeat= (currSubFeat<<2) | fNewBase;
				rNewWorkFeat= (currSubFeat<<2) | rNewBase;
				fNextSubFeat= (int) (fNewWorkFeat & subFeatureMask);
				rNextSubFeat= (int) (rNewWorkFeat & subFeatureMask);
				fShapeIndex	= (int) (fNewWorkFeat & 1023);
				rShapeIndex = (int) (rNewWorkFeat & 1023);
				fR0Weight	= 1;
				rR0Weight	= 1;
				if (index < l+flankLength+r0k-1) {
					fR0Weight = r0Alphas[(int) (fNewWorkFeat & r0FeatureMask)];
					rR0Weight = r0Alphas[(int) reverseComplement((rNewWorkFeat & r0FeatureMask), r0k)];
				}

				//Not in an active window.
				for (int currOffset=0; currOffset<startOffset; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset	= (index-currOffset-2)*nShapeClasses;
					fWeight			= fR0Weight;
					rWeight			= rR0Weight;
					//Adjust weights if shape features exist
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int csc=0; csc<nShapeClasses; csc++) {
							fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
							rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					//Update all positions
					for (int cgf=0; cgf<totFeatures; cgf++) {
						if (cgf<totNucFeatures) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totDinucFeatures) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totShapeFeatures) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
						}
						for (int cgf2=0; cgf2<=cgf; cgf2++) {
							fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
							rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
						}
					}
					fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
					revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
					fsWork[currOffset][fNextSubFeat] += fwdUpdate;
					rsWork[currOffset][rNextSubFeat] += revUpdate;
					
					//If shape features exist, add appropriate gradient components
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int csc=0; csc<nShapeClasses; csc++) {
							fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
									fwdUpdate*shapeFeatures[fShapeIndex][csc];
							rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
									revUpdate*shapeFeatures[rShapeIndex][csc];	
						}
						shapePosOffset += 4*k + 16*(k-1);
						for (int csc1=0; csc1<nShapeClasses; csc1++) {
							for (int csc2=0; csc2<=csc1; csc2++) {
								fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
								rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
							}
						}
						
						for (int csc=0; csc<nShapeClasses; csc++) {
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[shapePosOffset + csc][cgf][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][cgf][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf + 4*k;
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
						}
					}
				}
				
				//Active window.
				for (int currOffset=startOffset; currOffset<=endOffset; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					position		= index-currOffset;
					fNucPosOffset	= position*4 + fNewBase;
					rNucPosOffset	= position*4 + rNewBase;
					fDinucPosOffset = (position-1)*16 + (int) (fNewWorkFeat & 15);
					rDinucPosOffset = (position-1)*16 + (int) (rNewWorkFeat & 15);
					shapePosOffset	= (position-2)*nShapeClasses;
					fWeight			= fR0Weight;
					rWeight			= rR0Weight;
					if (position >= 1) {
						fWeight *= dinucAlphas[fDinucPosOffset];		
						rWeight *= dinucAlphas[rDinucPosOffset];			
					}
					fWeight *= nucAlphas[fNucPosOffset];
					rWeight *= nucAlphas[rNucPosOffset];
					//Adjust weights if shape features exist
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int csc=0; csc<nShapeClasses; csc++) {
							fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
							rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					//Update all positions
					for (int cgf=0; cgf<totFeatures; cgf++) {
						if (cgf<totNucFeatures) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totDinucFeatures) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totShapeFeatures) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
						}
						for (int cgf2=0; cgf2<=cgf; cgf2++) {
							fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
							rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
						}
					}
					fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
					revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
					fsWork[currOffset][fNextSubFeat]										+= fwdUpdate;
					rsWork[currOffset][rNextSubFeat]										+= revUpdate;						
					//Add to Nucleotide Gradients
					fsNucGradientsWork[fNucPosOffset][currOffset][fNextSubFeat]				+= fwdUpdate;						
					rsNucGradientsWork[rNucPosOffset][currOffset][rNextSubFeat]				+= revUpdate;
					//Add to Nucleotide Hessians
					fsHessianWork[fNucPosOffset][fNucPosOffset][currOffset][fNextSubFeat]	+= fwdUpdate;
					rsHessianWork[rNucPosOffset][rNucPosOffset][currOffset][rNextSubFeat]	+= revUpdate;

					//Do Shape features exist?
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						//Add to Shape Gradients
						for (int csc=0; csc<nShapeClasses; csc++) {
							fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
									fsOffsets[currOffset][currSubFeat]*fWeight*shapeFeatures[fShapeIndex][csc];
							rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
									rsOffsets[currOffset][currSubFeat]*rWeight*shapeFeatures[rShapeIndex][csc];
						}
						//Add to Shape Hessians
						shapePosOffset += 4*k + 16*(k-1);
						for (int csc1=0; csc1<nShapeClasses; csc1++) {
							for (int csc2=0; csc2<=csc1; csc2++) {
								fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
								rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
							}
							//Nuc/Shape Hessian Cross Terms
							fsHessianWork[shapePosOffset + csc1][fNucPosOffset][currOffset][fNextSubFeat] +=
									fwdUpdate*shapeFeatures[fShapeIndex][csc1];
							rsHessianWork[shapePosOffset + csc1][rNucPosOffset][currOffset][rNextSubFeat] +=
									revUpdate*shapeFeatures[rShapeIndex][csc1];
						}
						
						//Do dinucleotides exist?
						if (position>=1) {
							//Add to Dinucleotide Gradients
							fsDinucGradientsWork[fDinucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
							rsDinucGradientsWork[rDinucPosOffset][currOffset][rNextSubFeat] += revUpdate;
							//Add to Dinucleotide Hessians
							fDinucPosOffset += 4*k;
							rDinucPosOffset += 4*k;
							fsHessianWork[fDinucPosOffset][fDinucPosOffset][currOffset][fNextSubFeat]	+= fwdUpdate;
							rsHessianWork[rDinucPosOffset][rDinucPosOffset][currOffset][rNextSubFeat]	+= revUpdate;
							//Add to Nucleotide/Dinucleotide Hessian Cross Terms
							fsHessianWork[fDinucPosOffset][fNucPosOffset][currOffset][fNextSubFeat]		+= fwdUpdate;
							rsHessianWork[rDinucPosOffset][rNucPosOffset][currOffset][rNextSubFeat]		+= revUpdate;
							//Add to Shape/Dinucleotide Hessian Cross Terms
							for (int csc=0; csc<nShapeClasses; csc++) {
								fsHessianWork[shapePosOffset + csc][fDinucPosOffset][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc];
								rsHessianWork[shapePosOffset + csc][rDinucPosOffset][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc];
							}
							//Carry over all values from previous gradient
							//Carry over nucleotide features from all previous features into the nucleotide bin with the current features
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[fNucPosOffset][cgf][currOffset][fNextSubFeat]				+= fwdUpdate;
								rsHessianWork[rNucPosOffset][cgf][currOffset][rNextSubFeat]				+= revUpdate;
								fsHessianWork[fDinucPosOffset][cgf][currOffset][fNextSubFeat]			+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][cgf][currOffset][rNextSubFeat]			+= revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][cgf][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][cgf][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}
							}
							//Same, but with dinucleotide features
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k;
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat]			+= fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
								fsHessianWork[fDinucPosOffset][tempIdx][currOffset][fNextSubFeat]		+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][tempIdx][currOffset][rNextSubFeat]		+= revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}
							}
							//Same, but with shape features
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat] 			+= fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
								fsHessianWork[tempIdx][fDinucPosOffset][currOffset][fNextSubFeat]		+= fwdUpdate;
								rsHessianWork[tempIdx][rDinucPosOffset][currOffset][rNextSubFeat]		+= revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}		
							}
						} else {	//No Dinucleotides exist.
							//Carry over all values from previous gradient
							//Carry over nucleotide features from all previous features into the nucleotide bin with the current features 
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[fNucPosOffset][cgf][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[rNucPosOffset][cgf][currOffset][rNextSubFeat] += revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][cgf][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][cgf][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}
							}
							//Same, but with dinucleotide features
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k;
								fsHessianWork[fNucPosOffset][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[rNucPosOffset][tempIdx][currOffset][rNextSubFeat] += revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat] +=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat] +=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}
							}
							//Same, but with shape features
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat]	+= revUpdate;
								for (int csc=0; csc<nShapeClasses; csc++) {
									fsHessianWork[shapePosOffset+csc][tempIdx][currOffset][fNextSubFeat]+=
											fwdUpdate*shapeFeatures[fShapeIndex][csc];
									rsHessianWork[shapePosOffset+csc][tempIdx][currOffset][rNextSubFeat]+=
											revUpdate*shapeFeatures[rShapeIndex][csc];
								}		
							}
						}
					} else {	//No Shape features exist.
						//Do dinucleotides exist?
						if (position>=1) {
							//Add to Dinucleotide Gradients
							fsDinucGradientsWork[fDinucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
							rsDinucGradientsWork[rDinucPosOffset][currOffset][rNextSubFeat] += revUpdate;
							//Add to Dinucleotide Hessians
							fDinucPosOffset += 4*k;
							rDinucPosOffset += 4*k;
							fsHessianWork[fDinucPosOffset][fDinucPosOffset][currOffset][fNextSubFeat]	+= fwdUpdate;
							rsHessianWork[rDinucPosOffset][rDinucPosOffset][currOffset][rNextSubFeat]	+= revUpdate;
							//Add to Nucleotide/Dinucleotide Hessian Cross Terms
							fsHessianWork[fDinucPosOffset][fNucPosOffset][currOffset][fNextSubFeat]		+= fwdUpdate;
							rsHessianWork[rDinucPosOffset][rNucPosOffset][currOffset][rNextSubFeat]		+= revUpdate;
							//Carry over all values from previous gradient
							//Carry over nucleotide features from all previous features into the nucleotide bin with the current features
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[fNucPosOffset][cgf][currOffset][fNextSubFeat]				+= fwdUpdate;
								rsHessianWork[rNucPosOffset][cgf][currOffset][rNextSubFeat]				+= revUpdate;
								fsHessianWork[fDinucPosOffset][cgf][currOffset][fNextSubFeat]			+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][cgf][currOffset][rNextSubFeat]			+= revUpdate;
							}
							//Same, but with dinucleotide features
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k;
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat]			+= fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat]			+= revUpdate;
								fsHessianWork[fDinucPosOffset][tempIdx][currOffset][fNextSubFeat]		+= fwdUpdate;
								rsHessianWork[rDinucPosOffset][tempIdx][currOffset][rNextSubFeat]		+= revUpdate;
							}
							//Same, but with shape features
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat] += revUpdate;
								fsHessianWork[tempIdx][fDinucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[tempIdx][rDinucPosOffset][currOffset][rNextSubFeat] += revUpdate;
							}
						} else {	//No Dinucleotides exist.
							//Carry over all values from previous gradient
							//Carry over nucleotide features from all previous features into the nucleotide bin with the current features 
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[fNucPosOffset][cgf][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[rNucPosOffset][cgf][currOffset][rNextSubFeat] += revUpdate;
							}
							//Same, but with dinucleotide features
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k;
								fsHessianWork[fNucPosOffset][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[rNucPosOffset][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[tempIdx][fNucPosOffset][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[tempIdx][rNucPosOffset][currOffset][rNextSubFeat] += revUpdate;
							}
						}
					}
				}
				
				//Not in an active window.
				for (int currOffset=endOffset+1; currOffset<maxFrames; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset	= (index-currOffset-2)*nShapeClasses;
					fWeight			= fR0Weight;
					rWeight			= rR0Weight;
					//Adjust weights if shape features exist
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int csc=0; csc<nShapeClasses; csc++) {
							fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
							rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					} 
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int csc=0; csc<nShapeClasses; csc++) {
							fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
							rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					}
					//Update all positions
					for (int cgf=0; cgf<totFeatures; cgf++) {
						if (cgf<totNucFeatures) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totDinucFeatures) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totShapeFeatures) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
						}
						for (int cgf2=0; cgf2<=cgf; cgf2++) {
							fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
							rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
						}
					}
					fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
					revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
					fsWork[currOffset][fNextSubFeat] += fwdUpdate;
					rsWork[currOffset][rNextSubFeat] += revUpdate;
					//If shape features exist, add appropriate gradient components
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int csc=0; csc<nShapeClasses; csc++) {
							fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
									fwdUpdate*shapeFeatures[fShapeIndex][csc];
							rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
									revUpdate*shapeFeatures[rShapeIndex][csc];	
						}
						shapePosOffset += 4*k + 16*(k-1);
						for (int csc1=0; csc1<nShapeClasses; csc1++) {
							for (int csc2=0; csc2<=csc1; csc2++) {
								fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
								rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
							}
						}
						
						for (int csc=0; csc<nShapeClasses; csc++) {
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[shapePosOffset + csc][cgf][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][cgf][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf + 4*k;
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
						}
					}
				}
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradients[k][i][j] = fsNucGradientsWork[k][i][j];
						rsNucGradients[k][i][j] = rsNucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradients[k][i][j] = fsDinucGradientsWork[k][i][j];
						rsDinucGradients[k][i][j] = rsDinucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradients[k][i][j] = fsShapeGradientsWork[k][i][j];
						rsShapeGradients[k][i][j] = rsShapeGradientsWork[k][i][j];
					}
				}
			}
			for (int l=0; l<totFeatures; l++) {
				for (int k=0; k<totFeatures; k++) {
					for (int i=0; i<maxFrames; i++) {
						for (int j=0; j<maxSubFeatures; j++) {
							fsHessian[l][k][i][j] = fsHessianWork[l][k][i][j];
							rsHessian[l][k][i][j] = rsHessianWork[l][k][i][j];
						}
					}
				}				
			}
		}
		
		//Loop over last 2 bases to complete shape features
		int maxIndex = Math.max(2, r0k-1-flankLength);
		for (int index=l+2*flankLength; index<l+2*flankLength+maxIndex; index++) {	
			shapeStartOffset	= (index-k-1 < 0) ? 0 : index-k-1;
			shapeEndOffset		= (index-2 > maxFrames-1) ? maxFrames-1 : index-2;
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsWork[i][j] = 0;
					rsWork[i][j] = 0;
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradientsWork[k][i][j] = 0;
						rsNucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradientsWork[k][i][j] = 0;
						rsDinucGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradientsWork[k][i][j] = 0;
						rsShapeGradientsWork[k][i][j] = 0;
					}
				}
			}
			for (int l=0; l<totFeatures; l++) {
				for (int k=0; k<totFeatures; k++) {
					for (int i=0; i<maxFrames; i++) {
						for (int j=0; j<maxSubFeatures; j++) {
							fsHessianWork[l][k][i][j] = 0;
							rsHessianWork[l][k][i][j] = 0;
						}
					}
				}				
			}
			
			//Find new base
			int fNewBase = (int) fsRFlank & 3;
			int rNewBase = (int) rsRFlank & 3;
			fsRFlank >>= 2;
			rsRFlank >>= 2;
			for (int currSubFeat=0; currSubFeat<maxSubFeatures; currSubFeat++) {
				fNewWorkFeat= (currSubFeat<<2) | fNewBase;			//What is the new feature?
				rNewWorkFeat= (currSubFeat<<2) | rNewBase;
				fNextSubFeat= (int) (fNewWorkFeat & subFeatureMask);
				rNextSubFeat= (int) (rNewWorkFeat & subFeatureMask);
				fShapeIndex = (int) (fNewWorkFeat & 1023);
				rShapeIndex = (int) (rNewWorkFeat & 1023);
				fR0Weight	= 1;
				rR0Weight	= 1;
				if (index < l+flankLength+r0k-1) {
					fR0Weight = r0Alphas[(int) (fNewWorkFeat & r0FeatureMask)];
					rR0Weight = r0Alphas[(int) reverseComplement((rNewWorkFeat & r0FeatureMask), r0k)];					
				}
				for (int currOffset=0; currOffset<maxFrames; currOffset++) {
					if (fsOffsets[currOffset][currSubFeat]==0 && rsOffsets[currOffset][currSubFeat]==0)	continue;
					shapePosOffset	= (index-currOffset-2)*nShapeClasses;
					fWeight			= fR0Weight;
					rWeight			= rR0Weight;
					//Adjust weights if shape features exist
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						fShapeWeight = 0;
						rShapeWeight = 0;
						for (int csc=0; csc<nShapeClasses; csc++) {
							fShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[fShapeIndex][csc];
							rShapeWeight += shapeBetas[shapePosOffset + csc]*shapeFeatures[rShapeIndex][csc];
						}
						fWeight *= Math.exp(fShapeWeight);
						rWeight *= Math.exp(rShapeWeight);
					}
					//Update all positions
					for (int cgf=0; cgf<totFeatures; cgf++) {
						if (cgf<totNucFeatures) {
							fsNucGradientsWork[cgf][currOffset][fNextSubFeat] += fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsNucGradientsWork[cgf][currOffset][rNextSubFeat] += rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totDinucFeatures) {
							fsDinucGradientsWork[cgf][currOffset][fNextSubFeat] += fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsDinucGradientsWork[cgf][currOffset][rNextSubFeat] += rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
						}
						if (cgf<totShapeFeatures) {
							fsShapeGradientsWork[cgf][currOffset][fNextSubFeat] += fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
							rsShapeGradientsWork[cgf][currOffset][rNextSubFeat] += rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;									
						}
						for (int cgf2=0; cgf2<=cgf; cgf2++) {
							fsHessianWork[cgf][cgf2][currOffset][fNextSubFeat] += fsHessian[cgf][cgf2][currOffset][currSubFeat]*fWeight;
							rsHessianWork[cgf][cgf2][currOffset][rNextSubFeat] += rsHessian[cgf][cgf2][currOffset][currSubFeat]*rWeight;
						}
					}
					fwdUpdate = fsOffsets[currOffset][currSubFeat]*fWeight;
					revUpdate = rsOffsets[currOffset][currSubFeat]*rWeight;
					fsWork[currOffset][fNextSubFeat] += fwdUpdate;
					rsWork[currOffset][rNextSubFeat] += revUpdate;
					
					//If shape features exist, add appropriate gradient components
					if (currOffset>=shapeStartOffset && currOffset <= shapeEndOffset) {
						for (int csc=0; csc<nShapeClasses; csc++) {
							fsShapeGradientsWork[shapePosOffset + csc][currOffset][fNextSubFeat] +=
									fwdUpdate*shapeFeatures[fShapeIndex][csc];
							rsShapeGradientsWork[shapePosOffset + csc][currOffset][rNextSubFeat] +=
									revUpdate*shapeFeatures[rShapeIndex][csc];	
						}
						shapePosOffset += 4*k + 16*(k-1);
						for (int csc1=0; csc1<nShapeClasses; csc1++) {
							for (int csc2=0; csc2<=csc1; csc2++) {
								fsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][fNextSubFeat] +=
										fwdUpdate*shapeFeatures[fShapeIndex][csc1]*shapeFeatures[fShapeIndex][csc2];
								rsHessianWork[shapePosOffset+csc1][shapePosOffset+csc2][currOffset][rNextSubFeat] +=
										revUpdate*shapeFeatures[rShapeIndex][csc1]*shapeFeatures[rShapeIndex][csc2];
							}
						}
						for (int csc=0; csc<nShapeClasses; csc++) {
							for (int cgf=0; cgf<totNucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsNucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsNucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								fsHessianWork[shapePosOffset + csc][cgf][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][cgf][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totDinucFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsDinucGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsDinucGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k;
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat]	+= fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
							for (int cgf=0; cgf<totShapeFeatures; cgf++) {
								fwdUpdate = shapeFeatures[fShapeIndex][csc]*fsShapeGradients[cgf][currOffset][currSubFeat]*fWeight;
								revUpdate = shapeFeatures[rShapeIndex][csc]*rsShapeGradients[cgf][currOffset][currSubFeat]*rWeight;
								if (fwdUpdate==0 && revUpdate==0)	continue;
								tempIdx = cgf+4*k+16*(k-1);
								fsHessianWork[shapePosOffset + csc][tempIdx][currOffset][fNextSubFeat] += fwdUpdate;
								rsHessianWork[shapePosOffset + csc][tempIdx][currOffset][rNextSubFeat] += revUpdate;
							}
						}
					}
				}
			}
			for (int i=0; i<maxFrames; i++) {
				for (int j=0; j<maxSubFeatures; j++) {
					fsOffsets[i][j] = fsWork[i][j];
					rsOffsets[i][j] = rsWork[i][j];
				}
			}
			for (int k=0; k<totNucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsNucGradients[k][i][j] = fsNucGradientsWork[k][i][j];
						rsNucGradients[k][i][j] = rsNucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totDinucFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsDinucGradients[k][i][j] = fsDinucGradientsWork[k][i][j];
						rsDinucGradients[k][i][j] = rsDinucGradientsWork[k][i][j];
					}
				}
			}
			for (int k=0; k<totShapeFeatures; k++) {
				for (int i=0; i<maxFrames; i++) {
					for (int j=0; j<maxSubFeatures; j++) {
						fsShapeGradients[k][i][j] = fsShapeGradientsWork[k][i][j];
						rsShapeGradients[k][i][j] = rsShapeGradientsWork[k][i][j];
					}
				}
			}
			for (int l=0; l<totFeatures; l++) {
				for (int k=0; k<totFeatures; k++) {
					for (int i=0; i<maxFrames; i++) {
						for (int j=0; j<maxSubFeatures; j++) {
							fsHessian[l][k][i][j] = fsHessianWork[l][k][i][j];
							rsHessian[l][k][i][j] = rsHessianWork[l][k][i][j];
						}
					}
				}				
			}
		}
				
		Z = 0;
		for (int i=0; i<maxFrames; i++) {
			Z += Array.sum(fsOffsets[i]) + Array.sum(rsOffsets[i]);
		}
		//Sum over frames to get total gradient
		nucGradients	= new double[totNucFeatures];
		dinucGradients	= new double[totDinucFeatures];
		shapeGradients	= new double[totShapeFeatures];
		hessian			= (isNSBinding) ? new double[totFeatures+1][totFeatures+1] : new double[totFeatures][totFeatures];
		for (int cgf=0; cgf<totFeatures; cgf++) {
			for (int currOffset=0; currOffset<maxFrames; currOffset++) {				
				if (cgf<totNucFeatures) {
					nucGradients[cgf] += Array.sum(fsNucGradients[cgf][currOffset]) 
							+ Array.sum(rsNucGradients[cgf][currOffset]);
				}
				if (cgf<totDinucFeatures) {
					dinucGradients[cgf] += Array.sum(fsDinucGradients[cgf][currOffset]) 
							+ Array.sum(rsDinucGradients[cgf][currOffset]);	
				}
				if (cgf<totShapeFeatures) {
					shapeGradients[cgf] += Array.sum(fsShapeGradients[cgf][currOffset])
							+ Array.sum(rsShapeGradients[cgf][currOffset]);
				}
				for (int cgf2=0; cgf2<=cgf; cgf2++) {
					hessian[cgf][cgf2] += Array.sum(fsHessian[cgf][cgf2][currOffset]) 
							+ Array.sum(rsHessian[cgf][cgf2][currOffset]);
				}
			}
		}
		//Symmetrize
		for (int cgf=0; cgf<totFeatures; cgf++) {
			for (int cgf2=cgf+1; cgf2<totFeatures; cgf2++) {
				hessian[cgf][cgf2] = hessian[cgf2][cgf];
			}
		}
		
		fsOffsets			= null;
		rsOffsets			= null;
		fsWork				= null;
		rsWork				= null;
		fsNucGradients		= null;
		rsNucGradients		= null;
		fsDinucGradients	= null;
		rsDinucGradients	= null;
		fsShapeGradients	= null;
		rsShapeGradients	= null;
		fsNucGradientsWork	= null;
		rsNucGradientsWork	= null;
		fsDinucGradientsWork= null;
		rsDinucGradientsWork= null;
		fsShapeGradientsWork= null;
		rsShapeGradientsWork= null;
		fsHessian			= null;
		rsHessian			= null;
		fsHessianWork		= null;
		rsHessianWork		= null;
		System.gc();
	}
	
	@Override
	public double getZ() {
		return Z;
	}

	@Override
	public double[] getNucGradients() {
		return Array.clone(nucGradients);
	}

	@Override
	public double[] getDinucGradients() {
		return Array.clone(dinucGradients);
	}

	@Override
	public double[] getShapeGradients() {
		return Array.clone(shapeGradients);
	}
	
	@Override
	public double[][] getHessian() {
		return matrixClone(hessian);
	}
}
