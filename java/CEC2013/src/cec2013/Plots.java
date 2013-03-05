/******************************************************************************
 * Version: 1.0
 * Last modified on: 4 March, 2013 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
package cec2013;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

public final class Plots {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
	    /**********************************************************************
	     * Function 1 
	     *********************************************************************/
		System.out.println("Preparing files for Function 1...");
		Func f1 = new Funcs.FiveUnevenPeakTrap( new BoundsFn.ConstantBoundsFn( 1, new ClosedInterval.Double(0.0, 30.0) ) );
		Plot1D(f1, "plots/f1.dat", 0.1);
		
	    /**********************************************************************
	     * Function 2 
	     *********************************************************************/
		System.out.println("Preparing files for Function 2...");
		Func f2 = new Funcs.EqualMaxima( new BoundsFn.ConstantBoundsFn( 1, new ClosedInterval.Double(0.0, 1.0) ) );
		Plot1D(f2, "plots/f2.dat", 0.001);
		
	    /**********************************************************************
	     * Function 3 
	     *********************************************************************/
		System.out.println("Preparing files for Function 3...");
		Func f3 = new Funcs.UnevenDecreasingMaxima( new BoundsFn.ConstantBoundsFn( 1, new ClosedInterval.Double(0.0, 1.0) ) );
		Plot1D(f3, "plots/f3.dat", 0.001);
		
	    /**********************************************************************
	     * Function 4 
	     *********************************************************************/
		System.out.println("Preparing files for Function 4...");
		Func f4 = new Funcs.Himmelblau( new BoundsFn.ConstantBoundsFn( 2, new ClosedInterval.Double( -6.0, 6.0 ) ) );
		double [] x = new double[2]; x[0] = x[1] = 0.05;
		Plot2D(f4, "plots/f4.dat", x);

	    /**********************************************************************
	     * Function 5 
	     *********************************************************************/
		System.out.println("Preparing files for Function 5...");
		Func f5 = new Funcs.SixHumpCamelBack( new BoundsFn.ExplicitBoundsFn( 
				new ClosedInterval.Double( -1.9, 1.9 ), 
				new ClosedInterval.Double( -1.1, 1.1 ) ) );
		x[0] = 0.019; x[1] = 0.011;
		Plot2D(f5, "plots/f5.dat", x);

	    /**********************************************************************
	     * Function 6 
	     *********************************************************************/
		System.out.println("Preparing files for Function 6...");
		Func f6 = new Funcs.Shubert( new BoundsFn.ConstantBoundsFn( 2, new ClosedInterval.Double( -10.0, 10.0 ) ) );
		Plot2D(f6, "plots/f6.dat", new double [] {0.1, 0.1});

	    /**********************************************************************
	     * Function 7 
	     *********************************************************************/
		System.out.println("Preparing files for Function 7...");
		Func f7 = new Funcs.Vincent( new BoundsFn.ConstantBoundsFn( 2, new ClosedInterval.Double( 0.25, 10.0 ) ) );
		Plot2D(f7, "plots/f7.dat", new double [] {0.05, 0.05});

		/**********************************************************************
	     * Function 8 
	     *********************************************************************/
		System.out.println("Preparing files for Function 8...");
		Func f8 = new Funcs.ModifiedRastriginAll( new BoundsFn.ConstantBoundsFn( 2, new ClosedInterval.Double( 0.0, 1.0 ) ) );
		Plot2D(f8, "plots/f8.dat", new double [] {0.01, 0.01});

	    /**********************************************************************
	     * Composition functions F9-F10, 2-D data files for CF1-CF4
	     *********************************************************************/
		System.out.println("Preparing files for Composition Function 1...");
		Func cf1 = new Funcs.CF1( new BoundsFn.ConstantBoundsFn( 2, new ClosedInterval.Double( -5, 5 ) ) );
		Plot2D(cf1, "plots/CF1.dat", new double [] {0.01, 0.01});
		
		System.out.println("Preparing files for Composition Function 2...");
		Func cf2 = new Funcs.CF2( new BoundsFn.ConstantBoundsFn( 2, new ClosedInterval.Double( -5, 5 ) ) );
		Plot2D(cf2, "plots/CF2.dat", new double [] {0.01, 0.01});

		System.out.println("Preparing files for Composition Function 3...");
		Func cf3 = new Funcs.CF3( new BoundsFn.ConstantBoundsFn( 2, new ClosedInterval.Double( -5, 5 ) ) );
		Plot2D(cf3, "plots/CF3.dat", new double [] {0.01, 0.01});

		System.out.println("Preparing files for Composition Function 4...");
		Func cf4 = new Funcs.CF4( new BoundsFn.ConstantBoundsFn( 2, new ClosedInterval.Double( -5, 5 ) ) );
		Plot2D(cf4, "plots/CF4.dat", new double [] {0.01, 0.01});

	}//main
	
	public static void Plot1D(Func pFunc, String filename, final double stepsize) {
		try {
			/* Open the file */
			FileOutputStream fos   = new FileOutputStream(filename);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			BufferedWriter bw      = new BufferedWriter(osw);

			try {
				double []x = new double[pFunc.getDimension()];

				for (x[0]=pFunc.getBounds().get(0).getLower(); 
						x[0]<=pFunc.getBounds().get(0).getUpper(); 
						x[0]=x[0]+stepsize) {
					bw.write(x[0]+"\t"+pFunc.evaluate(x)+"\n");
				}
			} finally {
				/* Close the file */
				bw.close();

			}

		}catch (IOException e) {
			e.printStackTrace();
		}  
	}

	public static void Plot2D(Func pFunc, String filename, final double []stepsize) {
		try {
			/* Open the file */
			FileOutputStream fos   = new FileOutputStream(filename);
			OutputStreamWriter osw = new OutputStreamWriter(fos);
			BufferedWriter bw      = new BufferedWriter(osw);

			try {
				double []x = new double[pFunc.getDimension()];

				for (x[0]=pFunc.getBounds().get(0).getLower(); 
						x[0]<=pFunc.getBounds().get(0).getUpper(); 
						x[0]=x[0]+stepsize[0]) {
					for (x[1]=pFunc.getBounds().get(1).getLower(); 
							x[1]<=pFunc.getBounds().get(1).getUpper(); 
							x[1]=x[1]+stepsize[1]) {
						//System.out.println("0: "+x[0] + " 1: "+x[1]+" f: "+ pFunc.evaluate(x));
						bw.write(pFunc.evaluate(x)+"\t");
					}
					bw.newLine();
				}
			} finally {
				/* Close the file */
				bw.close();

			}

		}catch (IOException e) {
			e.printStackTrace();
		}       
	}
}
