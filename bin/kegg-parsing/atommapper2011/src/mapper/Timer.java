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

public class Timer
{
	private long startTime = 0;
	private long stopTime = 0;

	public Timer()
	{
	}

	public Timer(boolean init)
	{
		if (init)
			this.startTime = System.currentTimeMillis();
	}

	public void start()
	{
		this.startTime = System.currentTimeMillis();
	}

	public long getTime()
	{
		if (this.stopTime != 0)
			return this.stopTime - this.startTime;
		if (this.startTime != 0)
			return System.currentTimeMillis() - this.startTime;
		else
			return -1;
	}

	public void stop()
	{
		this.stopTime = System.currentTimeMillis();
	}

	public void reset()
	{
		this.startTime = 0;
		this.stopTime = 0;
	}
}

