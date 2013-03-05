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

public abstract class ClosedInterval
implements Comparable< ClosedInterval >, Cloneable
{
    public ClosedInterval() {}
    
    public static boolean isClosedInterval( double lower, double upper )
    {
	    return lower <= upper;
    }

    ///////////////////////////////
        
    public boolean isEmpty()
    {
        return getLength() == 0;
    }

    ///////////////////////////////
    
    public abstract double getLower();
    public abstract double getUpper();

    ///////////////////////////////
            
    public boolean contains( double x )
    {
	    return getLower() <= x && x <= getUpper();
    }

    public boolean overlaps( ClosedInterval other )
    {
	    return contains( other.getLower() ) || contains( other.getUpper() );
    }

    public boolean contains( ClosedInterval other )
    {
	    return contains( other.getLower() ) && contains( other.getUpper() );
    }

    public double getLength() 
    {
        double upper = getUpper();
        double lower = getLower();        
        return upper > lower ? upper - lower : 0;
    }

    public int hashCode( Object o ) {
    	return ( new java.lang.Double( getLower() ).hashCode() >> 13) ^ new java.lang.Double( getLower() ).hashCode();
    }
	
	public boolean equals( Object o )
    {
        return o instanceof ClosedInterval
                && getLower() == ((ClosedInterval)o).getLower()
                    && getUpper() == ((ClosedInterval)o).getUpper();
    }
    
    public String toString()
    {
        return "[ " + getLower() + ", " + getUpper() + " ]";
    }
            
    public int compareTo( ClosedInterval rhs )
    {
        if( getLower() < rhs.getLower() )
            return -1;
        else if( getLower() == rhs.getLower() )
        {
            if( getUpper() < rhs.getUpper() )
                return -1;
            else if( getUpper() == rhs.getUpper() )
                return 0;
            else
            	return 1;
        }
        else
        	return 1;        
    }    
    
    public Object clone()
    {
        return this;    
    }
    
    ///////////////////////////////
    
    public ClosedInterval.Double 
    intersection( ClosedInterval rhs )
    {
        ClosedInterval.Double result = new ClosedInterval.Double();

        if( overlaps( rhs ) )
        {
            double lower = getLower() > rhs.getLower() ? getLower() : rhs.getLower();
            double upper = getUpper() < rhs.getUpper() ? getUpper() : rhs.getUpper();
            result = new ClosedInterval.Double( lower, upper );
        }    

        return result;
    }
    
    ///////////////////////////////    

    public static class Double
    extends ClosedInterval
    {
        public double lower = 0;
        public double upper = -1;
        
        ///////////////////////////////
            
        public Double() {}

        public Double( double fl, double cl )
        {
        	if( !isClosedInterval( fl, cl ) )
        		throw new IllegalArgumentException( "closed interval expected, found [" + fl + ',' + cl + ']' );
        	
            lower = fl;
            upper = cl;
        }

        public Double( Double rhs )
        {
            lower = rhs.lower;
            upper = rhs.upper;
        }

        public double getLower() { return lower; }
        public double getUpper() { return upper; }

        public boolean equals( Object o )
        {
        	if( !( o instanceof  ClosedInterval.Double ) )
        		return false;
        	
        	ClosedInterval.Double rhs = (ClosedInterval.Double)o;
        	return lower == rhs.lower && upper == rhs.upper;
        }
    }
    
    ///////////////////////////////    
    
    public static class Int
    extends ClosedInterval
    {
        public int lower = 0;
        public int  upper = -1;
        
        ///////////////////////////////
            
        public Int() {}

        public Int( int fl, int cl )
        {
            lower = fl;
            upper = cl;
        }

        public Int( Int rhs )
        {
            lower = rhs.lower;
            upper = rhs.upper;
        }

        public double getLower() { return lower; }
        public double getUpper() { return upper; }

        public boolean equals( Object o )
        {
        	if( !( o instanceof  ClosedInterval.Long ) )
        		return false;
        	
        	ClosedInterval.Long rhs = (ClosedInterval.Long)o;
        	return lower == rhs.lower && upper == rhs.upper;
        }
    }
    
    ///////////////////////////////    

    public static class Long
    extends ClosedInterval
    {
        public long lower = 0;
        public long upper = -1;
        
        ///////////////////////////////
            
        public Long() {}

        public Long( long fl, long cl )
        {
            lower = fl;
            upper = cl;
        }

        public Long( Long rhs )
        {
            lower = rhs.lower;
            upper = rhs.upper;
        }

        public double getLower() { return lower; }
        public double getUpper() { return upper; }

        public boolean equals( Object o )
        {
        	if( !( o instanceof  ClosedInterval.Int ) )
        		return false;
        	
        	ClosedInterval.Int rhs = (ClosedInterval.Int)o;
        	return lower == rhs.lower && upper == rhs.upper;
        }
    }
    
    ///////////////////////////////
}

// End ///////////////////////////////////////////////////////////////
