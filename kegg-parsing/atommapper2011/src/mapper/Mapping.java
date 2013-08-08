//
//    Copyright 2011 Markus Heinonen 
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
package mapper;

import java.util.*;


// interface for mappings:
// this is implemented by HashMapping and AstarMapping
public interface Mapping
{
	public int size();
	public void extend(Atom lhs, Atom rhs);
	
	public Collection<Atom> getDomain();
	public Collection<Atom> getRange();
	public Mapping clone();
	
	public Atom getImage(Atom lhs);
	
	public String printMapping();
	
	public void setGCost(int c);
	public void setHCost(int c);
	public void setFeatDiff(double c);
	public int getFCost();
	public int getHCost();
	public int getGCost();
	public int getCost();
	public double getFeatDiff();
	public int costFunction();
	public int costFunction(Atom lhs, Atom rhs);
	public void remove(Atom lhs);
}
