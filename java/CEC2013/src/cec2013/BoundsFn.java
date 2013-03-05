/******************************************************************************
 * Version: 1.0
 * Last modified on: 4 March, 2013 
 *           
 * ACKNOWLEDGE: This class has been kindly provided by Dr. Jerry Swan.
 * 				We really want to thank Dr. Jerry Swan for his fruitful 
 * 				comments and suggestions that helped to improve, the design 
 * 				and the implementation of the competition in the Java 
 * 				programming language.  
 * ***************************************************************************/
package cec2013;

public interface BoundsFn {

	public int getDimension();
	public ClosedInterval.Double getBound( int dim );
	
	///////////////////////////////
	
	public static class ConstantBoundsFn implements BoundsFn
	{
		private int dim;
		private ClosedInterval.Double value;
		
		public ConstantBoundsFn( int dim, ClosedInterval.Double value ) {
			this.dim = dim;
			this.value = value;			
		}
		
		@Override
		public int getDimension() {
			return dim;
		}

		@Override
		public ClosedInterval.Double getBound(int dim) {
			return value;
		}
		
	}

	///////////////////////////////

	public static class ExplicitBoundsFn implements BoundsFn
	{
		ClosedInterval.Double [] bounds;
		
		public ExplicitBoundsFn( ClosedInterval.Double ... bounds ) {
			this.bounds = bounds;
		}
		
		@Override
		public int getDimension() {
			return bounds.length;
		}

		@Override
		public ClosedInterval.Double getBound(int dim) {
			return bounds[ dim ];
		}
	}
}

// End ///////////////////////////////////////////////////////////////
