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

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Stack;


/*
 * class to represent priorityqeueu as arrays indexed by cost
 * 
 */

public class PriorityArray //implements Iterable<AstarMapping>
{
//	class PAIter implements Iterator<AstarMapping>
//	{
//		Iterator<AstarMapping> it;
//		int ptr;
//		
//		public PAIter()
//		{
//			ptr = 0;
//			it = queue.get(ptr).iterator();
//		}
//		
//		public boolean hasNext()
//		{
//			if (it.hasNext() || ptr < ub-lb)
//				return true;
//			return false;
//		}
//		
//		public AstarMapping next()
//		{
//			if (it.hasNext())
//				return it.next();
//			
//			while (ptr < ub-lb) // still lists to go through
//			{
//				ptr++;
//				it = queue.get(ptr).iterator();
//				if (it.hasNext())
//					break;
//			}
//			
//			if (it.hasNext())
//				return it.next();
//			
//			return null;
//		}
//		
//		public void remove() {} // unsupported
//	}
//	
	
	
	private Stack<AstarMapping> completes;
	private ArrayList<LinkedList<AstarMapping>> queue;
	private int ub;
	private int lb;
	private int size;
	
	
	public PriorityArray(int ub)
	{
		this.lb = 0;
		this.ub = ub;
		this.size = 0;
		
		completes = new Stack<AstarMapping>();
		queue = new ArrayList<LinkedList<AstarMapping>>();
		
		// rest of the arrays for mappings with cost 'k'
		for (int i = 0; i <= ub; i++)
		{
			queue.add( new LinkedList<AstarMapping>() );
		}
	}
	
	public boolean isEmpty()
	{
		return size == 0;
	}
	
	public int size()
	{
		return size;
	}
	
	public void offer(AstarMapping m)
	{
		int c = m.getFCost();
		
		if (c < lb)
			lb = c;
		
		size++;
		
		if (m.complete())
			completes.add(m);
		else
		{
			// if better featscore, put in front, otherwise in back
			if (queue.get(c).isEmpty())
				queue.get(c).offer(m);
			else if (m.getFeatDiff() <= queue.get(c).peekFirst().getFeatDiff())
				queue.get(c).offerFirst(m);
			else
				queue.get(c).offerLast(m);
		}
	}
	
	public AstarMapping poll()
	{
		// check first whether complete mappings is not empty
		if (!completes.isEmpty())
		{
			size--;
			return completes.pop();
		}
		
		// otherwise go through rest of the lists
		for (int i = lb; i <= ub; i++)
		{
			if (!queue.get(i).isEmpty())
			{	
				size--;
				return queue.get(i).poll();
			}
			else
				lb++;
		}
		
		return null;
	}
	
	// remove everything with score ub or worse
	public void removeByScore(int c)
	{
		for (int i = ub; i >= c; i--)
		{
			size -= queue.get(i).size();
			queue.remove(i);
		}
		ub = c-1;
		
		System.gc();
	}

	// check that there are no maps below lb limit
	// sort by featdiff
	public void refresh()
	{
		for (int i = 0; i < lb; i++)
		{
			if (!queue.get(i).isEmpty())
			{
				System.out.println("oh noes, lb miss");
				lb = i;
				break;
			}
		}
		
		// only sort first list
		Collections.sort(queue.get(lb), new FeatDiffComparator());
		
//		for (int i = lb; i <= ub; i++)
//		{
//			Collections.sort(queue.get(i), new FeatDiffComparator());
//		}
	}
	
//	public Iterator<AstarMapping> iterator()
//	{
//		return new PAIter();
//	}
	
	public String toString()
	{
		String s = size + ":" + completes.size() + ":";
		for (int i = 0; i <= ub; i++)
			s += queue.get(i).size() + ",";
		return s;
	}
}
