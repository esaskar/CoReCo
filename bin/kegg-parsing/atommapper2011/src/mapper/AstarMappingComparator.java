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

import java.util.Comparator;

public class AstarMappingComparator implements Comparator<AstarMapping>
{
	public int compare(AstarMapping a, AstarMapping b)
	{
//		// (1a) complete mappings with smaller fcost are preferred
//		if (a.UMPA.length == 0 && b.UMPA.length == 0)
		if (a.size() == Astar.reacsize && b.size() == Astar.reacsize)
		{
			if (a.getFCost() < b.getFCost())
				return -1;
			else if (a.getFCost() > b.getFCost())
				return 1;

			return 0;
		}
		
		// (1b) complete mapping wins partial mapping
		else if (a.size() == Astar.reacsize )
			return -1;
		else if (b.size() == Astar.reacsize )
			return 1;
		
		// (2) cheaper fcost wins
		else if (a.getFCost() < b.getFCost())
			return -1;
		else if (a.getFCost() > b.getFCost())
			return 1;

		// (3) larger size mappings are better
		if (a.size() > b.size())
			return -1;
		else if (a.size() < b.size())
			return 1;
			
//		// (4) smaller gcost is preferred against large hcost
//		else if (a.map.getGCost() < b.map.getGCost())
//			return -1;
//		else if (a.map.getGCost() > b.map.getGCost())
//			return 1;
	
//		// (4) smaller gcost is preferred against large gcost
//		else if (a.map.getGCost() < b.map.getGCost())
//			return -1;
//		else if (a.map.getGCost() > b.map.getGCost())
//			return 1;

		// (4) smaller gcost is preferred against large hcost
		else if (a.getFeatDiff() < b.getFeatDiff())
			return -1;
		else if (a.getFeatDiff() > b.getFeatDiff())
			return 1;

//		// (4) bigger featdiff is preferred
//		else if (a.map.getFeatDiff() < b.map.getFeatDiff())
//			return -1;
//		else if (a.map.getFeatDiff() > b.map.getFeatDiff())
//			return 1;

//		// (3) larger size mappings are better
//		else if (a.map.size() > b.map.size())
//			return -1;
//		else if (a.map.size() < b.map.size())
//			return 1;
		
		return 0;
	}
}

class FeatDiffComparator implements Comparator<AstarMapping>
{
	public int compare(AstarMapping a, AstarMapping b)
	{
		if (a.getFeatDiff() < b.getFeatDiff())
			return -1;
		else if (a.getFeatDiff() > b.getFeatDiff())
			return 1;
		return 0;
	}
}

