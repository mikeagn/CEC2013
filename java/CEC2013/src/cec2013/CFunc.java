/******************************************************************************
 * Version: 1.0
 * Last modified on: 4 March, 2013 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
package cec2013;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.List;

public abstract class CFunc extends Func {

	protected List< Func > funcs_;
	//protected int dimension_;
	protected int nofunc_;
	protected double C_;// = 2000.0;
	protected double [] lambda_;
	protected double [] sigma_;
	protected double [] bias_;
	protected double [][] O_;
	protected double [][][] M_;
	protected double [] weight_;
	protected double [] fi_;
	protected double [] z_;
	protected double f_bias_;
	protected double [] fmaxi_;
	protected double [] tmpx_;

	public CFunc(BoundsFn boundsFn, final int nofunc  ) { 
		super(boundsFn);
		this.nofunc_ 		= nofunc;

		if (this.getDimension() < 1 || nofunc_ < 1) 
			throw new IllegalArgumentException();

		initializeCFunc();
	}

	private void initializeCFunc() {
		C_ 		= 2000.0;
		f_bias_ = 0;
		lambda_ = new double[nofunc_];
		sigma_  = new double[nofunc_];
		bias_   = new double[nofunc_];
		O_  	= new double[nofunc_][this.getDimension()];
		M_		= new double[nofunc_][this.getDimension()][this.getDimension()];

		weight_ = new double[nofunc_];
		fi_ 	= new double[nofunc_];
		z_  	= new double[this.getDimension()];
		fmaxi_  = new double[nofunc_];
		tmpx_   = new double[this.getDimension()];
	}
	
	void loadOptima(final String filename) 
	{
		File file = new File(filename); 
		try {
			LineNumberReader reader = new LineNumberReader( new FileReader( file ) );
			try {
				String buffer;	
				try {
					for (int i=0; i<this.nofunc_; ++i) {
						buffer = reader.readLine();
						
						double tmp = 0;
						String[] numbers = buffer.split("\t");
						for( int j = 0; j < this.getDimension(); ++j ) {
							tmp = Double.parseDouble(numbers[j].trim());
							O_[i][j] = tmp;
						}
					}
				} catch (NumberFormatException e) {
					System.err.println("Error: loadOptima, NumberFormatException: " + e.toString());
					e.printStackTrace();
				} catch (IOException e) {
					System.err.println("Error: loadOptima, IOException: " + e.toString());
					e.printStackTrace();
				}
			} finally {
				try {
					reader.close();
				} catch (IOException e) {
					System.err.println("Error: loadOptima, Can not close file: IOException: " + e.toString());
					e.printStackTrace();
				}			
			}
		} catch (FileNotFoundException e) {
			System.err.println("Error: loadOptima, Can not find file: " + filename +e.toString());
			e.printStackTrace();
		} 
	}
	
	
	void loadRotationMatrix(final String filename)
	{	
		File file = new File(filename); 
		try {
			LineNumberReader reader = new LineNumberReader( new FileReader( file ) );
			try {
				String buffer;	
				try {
				    double tmp = -1;
				    for (int i=0; i<nofunc_; ++i) {
//				      System.out.println("Matrix: "+i);
				        for (int j=0; j<this.getDimension(); ++j) {
							buffer = reader.readLine();
							String[] numbers = buffer.split("\t");
				            for (int k=0; k<this.getDimension(); ++k) {
				            	tmp = Double.parseDouble(numbers[k].trim());
				                M_[i][j][k] = tmp;
//				                System.out.print(M_[i][j][k]+"\t");
				            }
//				            System.out.println("");
				        }
//				        System.out.println("");
				    }
				} catch (NumberFormatException e) {
					System.err.println("Error: loadRotationMatrix, NumberFormatException: " + e.toString());
					e.printStackTrace();
				} catch (IOException e) {
					System.err.println("Error: loadRotationMatrix, IOException: " + e.toString());
					e.printStackTrace();
				}
			} finally {
				try {
					reader.close();
				} catch (IOException e) {
					System.err.println("Error: loadRotationMatrix, Can not close file: IOException: " + e.toString());
					e.printStackTrace();
				}			
			}
		} catch (FileNotFoundException e) {
			System.err.println("Error: loadRotationMatrix, Can not find file: " + filename +e.toString());
			e.printStackTrace();
		} 
	}

	void calculateWeights(final double [] x)
	{
		double sum = 0, maxi = (-Double.POSITIVE_INFINITY), maxindex = 0;

		for (int i=0; i<nofunc_; ++i) {
			sum = 0.0;
			for (int j=0; j<this.getDimension(); ++j) {
				sum += ( x[j] - O_[i][j] ) * ( x[j] - O_[i][j] );
			}
			weight_[i] = Math.exp( -sum/(2.0 * this.getDimension() * sigma_[i] * sigma_[i]) );
			if (i==0) { maxi = weight_[i]; }
			if (weight_[i] > maxi) {
				maxi = weight_[i];
				maxindex = i;
			}
			//maxi = max(maxi, weight_[i]);
		}
		sum = 0.0;
		for (int i=0; i<nofunc_; ++i) {
			//if (weight_[i] != maxi) {
			if (i != maxindex) {
				weight_[i] *= (1.0 - Math.pow(maxi, 10.0));
			}
			sum += weight_[i];
		}
		for (int i=0; i<nofunc_; ++i) {
			if (sum == 0.0) {
				weight_[i] = 1.0/(double)nofunc_;
			} else {
				weight_[i] /= sum;
			}
		}
//		for (int i=0; i<nofunc_; ++i) {
//			System.out.print(weight_[i] +"\t");
//		}
//		System.out.println("");
	}

	void initRotmatIdentity()
	{
		for (int i=0; i<nofunc_; ++i) {
//			System.out.println("Matrix: "+ i);
			for (int j=0; j<this.getDimension(); ++j) {
				for (int k=0; k<this.getDimension(); ++k) {
					M_[i][j][k] = (j==k ? 1 : 0 );
//					System.out.print(M_[i][j][k] + "\t");
				}
//				System.out.println("");
			}
//			System.out.println("");
		}
	}
	void initOptimaRandomly()
	{
		List< ClosedInterval.Double > bounds = this.getBounds();
		for (int i=0; i< nofunc_; ++i) {
			for (int j=0; j< this.getDimension(); ++j) {
				//O_[i][j] = lbound_[j] + (ubound_[j] - lbound_[j]) * Math.random();//nextUniform(0.0,1.0);
				O_[i][j] = bounds.get(j).getLower() + 
						(bounds.get(j).getUpper() - bounds.get(j).getLower()) * Math.random();
				//	          System.out.print(O_[i][j] +"\t");
			}
			//	      System.out.println("");
		}
	}
	void transformToZ(final double []x, final int index)
	{
		/* Calculate z_i = (x - o_i)/\lambda_i */
		for (int i=0; i<this.getDimension(); ++i) {
			tmpx_[i] = (x[i] - O_[index][i])/lambda_[index];
		}
		/* Multiply z_i * M_i */
		for (int i=0; i<this.getDimension(); ++i) {
			z_[i] = 0;
			for (int j=0; j<this.getDimension(); ++j) {
				/* in MATLAB: M.M1*tmpx' */
				//z_[i] += M_[index][i][j] * tmpx_[j];

				/* in MATLAB: tmpx*M.M1 */
				z_[i] += M_[index][j][i] * tmpx_[j];
			}
			//System.out.println("i: "+i+" "+tmpx_[i]+" "+x[i]+" "+O_[index][i]+" "+lambda_[index]+" "+z_[i]);
		}
	}

	void transformToZNoshift(final double []x, final int index)
	{
		/* Calculate z_i = (x - o_i)/\lambda_i */
		for (int i=0; i<this.getDimension(); ++i) {
			//tmpx_[i] = (x[i] - O_[index][i])/lambda_[index];
			tmpx_[i] = (x[i])/lambda_[index];
		}
		/* Multiply z_i * M_i */
		for (int i=0; i<this.getDimension(); ++i) {
			z_[i] = 0;
			for (int j=0; j<this.getDimension(); ++j) {
				/* in MATLAB: M.M1*tmpx' */
				//z_[i] += M_[index][i][j] * tmpx_[j];

				/* in MATLAB: tmpx*M.M1 */
				z_[i] += M_[index][j][i] * tmpx_[j];
			}
		}
	}
	
	void CalculateFMaxi()
	{
		/* functions */
		double [] x5 = new double[this.getDimension()];
		for (int i=0; i<this.getDimension(); ++i) { x5[i] = 5 ; }
		int index=0;
		for( Func f : funcs_ ) {
			transformToZNoshift(x5, index);
//			for (int i=0; i<z_.length; i++)
//				System.out.println("z_: "+this.z_[i]);
			fmaxi_[index] = f.cfuncEvaluate(this.z_);
//			System.out.println("FMAXI: "+index+" : "+fmaxi_[index]);
			index++;
		}
	}
	
	double evaluateInner_(final double []x)
	{
	    double result = 0;
	    calculateWeights(x);
	    int index=0;
	    for( Func f : funcs_ ) {
	        transformToZ(x, index);
	        fi_[index] = f.cfuncEvaluate(this.z_);
//	        System.out.println("Func: "+i+" : "+fi_[i]);
	        index++;
	    }
	    for (int i=0; i<nofunc_; ++i) {
	        result += weight_[i]*( C_ * fi_[i] / fmaxi_[i] + bias_[i] );
	    }
	    //Assuming maximisation
	    return -1.0*result + f_bias_;
	}

	public abstract double doEvaluate(double [] x);

}

// End ///////////////////////////////////////////////////////////////
