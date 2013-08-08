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
import java.io.*;

public class Molecule implements Comparable<Molecule>
{
	private String id;
	private int molnum;
	private List<Atom> atoms = new ArrayList<Atom>();
	private List<Bond> bonds = new ArrayList<Bond>();
	private Atom maxDistAtom = null;
	
	
	public Atom getMaxDistAtom()
	{
		return maxDistAtom;
	}

	public void setMaxDistAtom(Atom maxDistAtom)
	{
		this.maxDistAtom = maxDistAtom;
	}

	public void computeMaxDistAtom()
	{
		// Find the most extremal atom (largest distance to any other node) for all molecules
		int maxdist = Integer.MIN_VALUE;
		Atom maxdistatom = null;
		for (Atom a : getAtoms())
		{
			int dist = traverse(a);
			
			if (dist > maxdist)
			{
				maxdist = dist;
				maxdistatom = a;
			}
		}
		
		this.maxDistAtom = maxdistatom;
	}
	
	private int traverse(Atom start)
	{
		int[] dist = new int[start.getParent().getAtoms().size()];
		LinkedList<Atom> Q = new LinkedList<Atom>();
		Q.addFirst(start);
		Atom node;
		dist[start.getMolecularId()] = 0;
		
		while (!Q.isEmpty())
		{
			node = Q.removeLast();
			
			for (Atom ne : node.getAtomNeighbors())
			{
				if (dist[ne.getMolecularId()] == 0)
				{
					Q.addFirst(ne);
					dist[ne.getMolecularId()] = dist[node.getMolecularId()] + 1;
				}
			}
		}
		
		int max = 0;
		for (int x : dist)
			if (x > max)
				max = x;
		
		return max;
	}
	
	public Molecule(int index, String id)
	{
		this.id = id;
		this.molnum = index;
	}

	public Molecule(Collection<Atom> at, Collection<Bond> bo, boolean setparent)
	{
		this.atoms = (ArrayList)at;
		this.bonds = (ArrayList)bo;
		if (setparent)
		{
			for (Atom a : at)
			{
				a.setParent(this);
			}
			for (Bond b : bo)
			{
				b.setParent(this);
			}
		}
	}

	public void addAtom(Atom a)
	{
		a.setParent(this);
		this.atoms.add(a);
	}

	public void addAtoms(Collection<Atom> at)
	{
		int i = this.atoms.size();
		this.atoms.addAll(at);
		for (Atom a : at)
		{
			a.setMolecularId(i++);
			a.setParent(this);
		}
	}

	public void addBond(Bond b)
	{
		b.setParent(this);
		this.bonds.add(b);
	}

	public void addBonds(Collection<Bond> bo)
	{
		this.bonds.addAll(bo);
		for (Bond b : bo)
		{
			b.setParent(this);
		}
	}

	public HashMap<String, Integer> atomSpectrum()
	{
		HashMap<String, Integer> h = new HashMap<String, Integer>();
		for (Atom a : this.atoms)
		{
			String s = a.getSymbol();
			if (h.containsKey(s))
			{
				int i = (int) h.get(s);
				i = i + 1;
				h.put(s, i);
			} else
			{
				h.put(s, 1);
			}
		}
		return h;
	}

	public HashMap<String, Integer> bondSpectrum()
	{
		HashMap<String, Integer> h = new HashMap<String, Integer>();
		for (Bond b : this.bonds)
		{
			String s = b.getBondString();
			
			if (h.containsKey(s))
				h.put(s, h.get(s)+1);
			else
				h.put(s, +1);
		}
		
		return h;
	}

	public Collection<Atom> getAtoms()
	{
		return this.atoms;
	}

	public void setMolnum(int m)
	{
		this.molnum = m;
	}

	public Bond getBond(Atom at1, Atom at2)
	{
		return at1.getBond(at2);
	}

	public Collection<Bond> getBonds()
	{
		return this.bonds;
	}

	public String getId()
	{
		return this.id;
	}

	public int getMolNum()
	{
		return this.molnum;
	}

	public int numAtoms()
	{
		return this.atoms.size();
	}

	public double atomSpectrumSimilarity(Molecule other)
	{
		int differing = 0;
		int total = 0;
		
		Map<String, Integer> ts = this.atomSpectrum();
		Map<String, Integer> os = other.atomSpectrum();
		
		Collection<String> keys = new HashSet<String>();
		keys.addAll(ts.keySet());
		keys.addAll(os.keySet());
		
		for (String key : keys)
		{
			int diff = 0;
			
			if (ts.containsKey(key))
			{
				diff += ts.get(key);
				total += ts.get(key);
			}
			if (os.containsKey(key))
			{
				diff -= os.get(key);
				total += os.get(key);
			}

			diff = Math.abs(diff);
			differing += diff;
		}
		
//		int total = ts.size() + os.size();
		
		// mtsvien lkm ja ei-mtsvien lkm
		// esim. 3 ja 4 atomia, 2 samoja, niin m=2, ei-m=3  -->  4/7?
		
		double result = 1.0 * (total-differing) / total;
		
		return result;
	}
	
	public double bondSpectrumSimilarity(Molecule other)
	{
		int differing = 0;
		int total = 0;
		
		Map<String, Integer> ts = this.bondSpectrum();
		Map<String, Integer> os = other.bondSpectrum();
		
		Collection<String> keys = new HashSet<String>();
		keys.addAll(ts.keySet());
		keys.addAll(os.keySet());
		
		for (String key : keys)
		{
			int diff = 0;
			
			if (ts.containsKey(key))
			{
				diff += ts.get(key);
				total += ts.get(key);
			}
			if (os.containsKey(key))
			{
				diff -= os.get(key);
				total += os.get(key);
			}

			diff = Math.abs(diff);
			differing += diff;
		}
		
//		int total = ts.size() + os.size();
		double result = 1.0 * (total-differing) / total;
		return result;
		
		// mtsvien lkm ja ei-mtsvien lkm
		// esim. 3 ja 4 atomia, 2 samoja, niin m=2, ei-m=3  -->  4/7?
	}
	
	public String toString()
	{
		StringBuffer sbf = new StringBuffer();
		sbf.append(this.numAtoms() + " atoms: ");
		for (Atom a : this.atoms)
			sbf.append(a.getSymbol());
		return sbf.toString();
	}

	public static Collection<Bond> verticesSquared(Collection<Atom> atoms)
	{
		Collection<Bond> atoms_2 = new HashSet<Bond>();
		for (Atom a : atoms)
		{
			for (Atom b : atoms)
			{
				atoms_2.add(new Bond(a, b, 1, false));
			}
		}
		return atoms_2;
	}

	public void read(String filename, String path, int reac_id_start) throws IOException
	{
		List<String> lines = new ArrayList<String>();
		
		Scanner sc = new Scanner(new File(path + filename));
		while (sc.hasNextLine())
			lines.add(sc.nextLine());
		
		int reac_id = reac_id_start; // start reaction's internal numberin from predefined number
		int id = 0; // start molecule's internal numbering from zero

		int atomcount = Integer.parseInt(lines.get(3).substring(0, 3).trim());
		int bondcount = Integer.parseInt(lines.get(3).substring(3, 6).trim());
		
		// contains information whether atom is hydrogen or not
		int[] validatoms = new int[atomcount];
		
		// atom block
		for (int i = 4; i < 4 + atomcount; i++)
		{
			String line = lines.get(i);
			String[] words = line.split("\\s+");
			String symbol = words[4].trim();
			
			// don't take hydrogens
			if (!symbol.equals("H") && !symbol.equals("H+"))
			{
				Atom a = new Atom(id++, symbol);
				a.setReactionalId(reac_id++);
				a.setParent(this);
				atoms.add(a);
				validatoms[i-4] = id-1;
			}
			else
				// all hydrogens into rawlist
				validatoms[i-4] = -1;
		}
		
		// bond block, ignore hydrogen bonds
		for (int i = 4 + atomcount; i < 4 + atomcount + bondcount; i++)
		{
			String line = lines.get(i);
			int source = Integer.parseInt(line.substring(0, 3).trim())-1; // numbering correction
			int target = Integer.parseInt(line.substring(3, 6).trim())-1; // mol-files start from 1, we from 0
			int type   = Integer.parseInt(line.substring(6, 9).trim());
			
			// hydrogen bonds are not created
			if (validatoms[source] == -1 || validatoms[target] == -1)
				continue;
			
			Bond b = new Bond(atoms.get(validatoms[source]), atoms.get(validatoms[target]), type, true);
			bonds.add(b);
		}
	}
	
	public int compareTo(Molecule other)
	{
		if (this.atoms.size() > other.getAtoms().size())
			return -1;
		else if (this.atoms.size() < other.getAtoms().size())
			return 1;

		return 0;
	}
	
}

