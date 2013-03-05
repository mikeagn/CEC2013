/******************************************************************************
 * Version: 1.0
 * Last modified on: 4 March, 2013 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
package cec2013;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ExampleUsage {

	public static void main( String [] args ) {
		/** Always Initialise a CEC2013 object **/
		CEC2013 comp = new CEC2013();
		/** Get the function list **/
		List< Func > funcs = comp.getFunctions();
		
		/*********************************************************************
		 * Example usage on how to iterate through the functions 
		 *********************************************************************/
		int index = 1;
		for( Func f : funcs ) {
			
			double []x = new double[f.getDimension()];
			Arrays.fill(x, 1);
			
			System.out.println( "f"+index+"(1..1)= " + f.evaluate(x) );
			index++;
		}
		/*********************************************************************
		 * Example usage on how to get iterate through the functions and 
		 * get the optimisation box (bounds) of each function 
		 *********************************************************************/
//		System.out.println("Optimisation bounds for each function: ");
//		List< ClosedInterval.Double > bounds;
//		index = 1;
//		for( Func f : funcs ) {
//			bounds = f.getBounds();
//			System.out.println("\nF"+index + ": ");
//			for (int i=0; i<f.getDimension(); ++i) {
//				System.out.println("i: "+i+" ["
//						+ bounds.get(i).getLower() + "," 
//						+ bounds.get(i).getUpper() + "]");
//			}
//			index++;
//		}
		
		/*********************************************************************
		 * Demonstration on using howManyGlobalOptimaInPopulation method
		 *********************************************************************/
		/* Create and initialise a population inside bounds */
		final int NP = 30; // population size
		ArrayList<ArrayList<Double>> population;
		List< ClosedInterval.Double > bounds;
		double accuracy = 0.001;

		/* Be careful findex has to start from value equal to 1 */
		int findex = 1;
		for( Func f : funcs ) {
			bounds = f.getBounds();
			population = new ArrayList<ArrayList<Double>>();

			/* Initialise the population randomly within bounds */
			for (int i=0; i<NP; ++i) {
				ArrayList<Double> indi = new ArrayList<Double>();
				for(int j=0; j<f.getDimension(); ++j) {
					Double tmp = bounds.get(j).getLower() + Math.random() * 
							(bounds.get(j).getUpper() - bounds.get(j).getLower());
					indi.add(tmp);
				}
				population.add(indi);
			}
			/* Count how many global optima exist in the population */
			ArrayList<ArrayList<Double>> seeds = new ArrayList<ArrayList<Double>>();
			int numberOfOptimaFound = CEC2013.howManyGlobalOptimaInPopulation(
							population, seeds, findex, accuracy, CEC2013.getRho(findex));
			
			System.out.println("In the current population there exist "+ 
					numberOfOptimaFound + " global optimizers.");
			
			if (numberOfOptimaFound != 0) {
				/* Print seeds */
				for(ArrayList<Double> indi : seeds){
					for (int i=0; i<indi.size(); ++i) 
						System.out.print(indi.get(i).doubleValue()+"\t");
					System.out.println(" Fitness: "+ f.evaluate(indi));
				}
			}
			findex++;
		}
	}
}
