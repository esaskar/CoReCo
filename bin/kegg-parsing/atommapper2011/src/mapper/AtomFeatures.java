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

public class AtomFeatures
{
	// atomin piirrevektori konteksteittain
	// wiener[2] = Map<atom,int>
	private Map<Integer,Integer> wiener;
	private Map<Integer,Map<Atom,Integer>> wienerdirs;
	private Map<Integer,String> ad;
	private Map<Integer,Map<Atom,String>> addirs;
	private Map<Integer,String> bd;
	private Map<Integer,Map<Atom,String>> bddirs;
	
	public AtomFeatures()
	{
		wiener = new HashMap<Integer,Integer>();
		wienerdirs = new HashMap<Integer,Map<Atom,Integer>>();
		ad = new HashMap<Integer,String>();
		addirs = new HashMap<Integer,Map<Atom,String>>();
		bd = new HashMap<Integer,String>();
		bddirs = new HashMap<Integer,Map<Atom,String>>();
	}
	// kahden atomin v�linen piirrevektori on konteksteittain eri suuntien samuus
	// atom 3 <-> atom 9': ( (1,1),(1,1),(0,0),(0,0) )
	//                   => (1,1,0,0)
	
	public void setWiener(int context, int val)
	{
		wiener.put(context, val);
	}
	
	public void setWienerDir(int context, Atom dir, int val)
	{
		if (!wienerdirs.containsKey(context))
			wienerdirs.put(context, new HashMap<Atom,Integer>());
		
		wienerdirs.get(context).put(dir, val);
	}

	public void setAd(int context, String val)
	{
		ad.put(context, val);
	}
	
	public void setAdDir(int context, Atom dir, String val)
	{
		if (!addirs.containsKey(context))
			addirs.put(context, new HashMap<Atom,String>());
		
		addirs.get(context).put(dir, val);
	}
	
	public void setBd(int context, String val)
	{
		bd.put(context, val);
	}
	
	public void setbdDir(int context, Atom dir, String val)
	{
		if (!bddirs.containsKey(context))
			bddirs.put(context, new HashMap<Atom,String>());
		
		bddirs.get(context).put(dir, val);
	}
}
