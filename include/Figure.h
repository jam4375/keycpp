// Figure.h -- Matlab-like interface to gnuplot through gnuplot-cpp

#ifndef FIGURE_H_
#define FIGURE_H_

#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <string>
#include <sstream>
#include "gnuplot_i.h"
#include "Matrix.h"

namespace keycpp
{
	class FigureException : public std::runtime_error
	{
		public:
			FigureException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	class Figure
	{
	private:
		Gnuplot g;
		bool hold_on_bool = false;
		int num_plots;
		matrix<double> colors;
		std::vector<std::vector<double> > x_plot_data;
		std::vector<std::vector<double> > y_plot_data;
		std::vector<std::string> plot_format;
		std::vector<double> plot_linewidth;
		std::vector<double> plot_markersize;
		std::vector<double> plot_val;
		std::vector<std::string> legend_entries;
		std::string m_xlabel;
		std::string m_ylabel;
		std::string m_title;
		bool grid_on_bool = false;
		int fontsize = 10;
	
	public:
		Figure();
		template<class U, class T> void plot(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1);
		template<class U, class T> void plot(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2);
		template<class U, class T> void plot(std::vector<U> x, std::vector<T> y, std::string arguments = "", double val = -1, double lw = 2, double ps = 1.5, std::string legend_entry = "");
		void xlabel(std::string xlabel_text);
		void ylabel(std::string ylabel_text);
		void title(std::string title_text);
		void grid_on();
		void grid_off();
		void hold_on();
		void hold_off();
		void legend(std::initializer_list<std::string> lst);
		void replot_all();
		void setFontsize(int p_fontsize) {fontsize = p_fontsize;};
		int getFontsize() {return fontsize;};
		void set(std::string property, double val);
	};
	
	inline Figure::Figure() try : g("lines"), num_plots(0)
	{
		colors = {{0.0,0.0,1.0},
			{1.0,0.0,0.0},
			{0.0,1.0,0.0},
			{0.0,0.0,0.172413793103448},
			{1.0,0.103448275862069,0.724137931034483},
			{1.0,0.827586206896552,0.0},
			{0.0,0.344827586206897,0.0},
			{0.517241379310345,0.517241379310345,1.0},
			{0.620689655172414,0.310344827586207,0.275862068965517},
			{0.0,1.0,0.758620689655172},
			{0.0,0.517241379310345,0.586206896551724},
			{0.0,0.0,0.482758620689655},
			{0.586206896551724,0.827586206896552,0.310344827586207},
			{0.965517241379310,0.620689655172414,0.862068965517241},
			{0.827586206896552,0.0689655172413793,1.0},
			{0.482758620689655,0.103448275862069,0.413793103448276},
			{0.965517241379310,0.0689655172413793,0.379310344827586},
			{1.0,0.758620689655172,0.517241379310345},
			{0.137931034482759,0.137931034482759,0.0344827586206897},
			{0.551724137931035,0.655172413793103,0.482758620689655}};
	}
	catch(GnuplotException ge)
	{
		std::cout << ge.what() << std::endl;
	}
	
	template<class U, class T> void Figure::plot(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1)
	{
		std::transform(format.begin(), format.end(), format.begin(), ::tolower);
		std::transform(property1.begin(), property1.end(), property1.begin(), ::tolower);
		
		double lw = 2;
		double ps = 1.5;
		if(!property1.empty())
		{
			if(property1.find("linewidth") != std::string::npos)
			{
				property1.erase(property1.find("linewidth"),9);
				if(!property1.empty())
				{
					throw FigureException("Unknown property string in Figure!");
				}
				lw = val1;
				plot(x,y,format,-1,lw,ps);
			}
			else if(property1.find("markersize") != std::string::npos)
			{
				property1.erase(property1.find("markersize"),10);
				if(!property1.empty())
				{
					throw FigureException("Unknown property string in Figure!");
				}
				ps = val1;
				plot(x,y,format,-1,lw,ps);
			}
		}
		else
		{
			throw FigureException("Invalid property while plotting in Figure!");
		}
	}
	
	template<class U, class T> void Figure::plot(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2)
	{
		std::transform(format.begin(), format.end(), format.begin(), ::tolower);
		std::transform(property1.begin(), property1.end(), property1.begin(), ::tolower);
		std::transform(property2.begin(), property2.end(), property2.begin(), ::tolower);
		double lw = 2;
		double ps = 1.5;
		bool lw_found = false, ps_found = false;
		if(!property1.empty())
		{
			if(property1.find("linewidth") != std::string::npos)
			{
				property1.erase(property1.find("linewidth"),9);
				if(!property1.empty())
				{
					throw FigureException("Unknown property string in Figure!");
				}
				lw = val1;
				lw_found = true;
			}
			else if(property1.find("markersize") != std::string::npos)
			{
				property1.erase(property1.find("markersize"),10);
				if(!property1.empty())
				{
					throw FigureException("Unknown property string in Figure!");
				}
				ps = val1;
				ps_found = true;
			}
		}
		else
		{
			throw FigureException("Invalid property while plotting in Figure!");
		}
		
		if(!property2.empty())
		{
			if(property2.find("linewidth") != std::string::npos)
			{
				property2.erase(property2.find("linewidth"),9);
				if(!property2.empty())
				{
					throw FigureException("Unknown property string in Figure!");
				}
				lw = val2;
				if(lw_found == true)
				{
					throw FigureException("Property is specified multiple times while plotting!");
				}
			}
			else if(property2.find("markersize") != std::string::npos)
			{
				property2.erase(property2.find("markersize"),10);
				if(!property2.empty())
				{
					throw FigureException("Unknown property string in Figure!");
				}
				ps = val2;
				if(ps_found == true)
				{
					throw FigureException("Property is specified multiple times while plotting!");
				}
			}
		}
		else
		{
			throw FigureException("Invalid property while plotting in Figure!");
		}
		plot(x,y,format,-1,lw,ps);
	}
	
	template<class U, class T> void Figure::plot(std::vector<U> x, std::vector<T> y, std::string arguments, double val, double lw, double ps, std::string legend_entry)
	{
		std::string format;
		try
		{
			std::transform(arguments.begin(), arguments.end(), arguments.begin(), ::tolower);
			format = arguments;
			num_plots++;
			int lt = -1;
			int pt = -1;
			std::string color_str = "";
			if(!arguments.empty())
			{
				if(arguments.find("linewidth") != std::string::npos)
				{
					arguments.erase(arguments.find("linewidth"),9);
					lw = val;
					lt = 1;
					if(!arguments.empty())
					{
						throw FigureException("Unknown property string in Figure!");
					}
				}
				else
				{
					if(arguments.find("--") != std::string::npos && lt < 0)
					{
						arguments.erase(arguments.find("--"),2);
						lt = 2;
					}
					else if(arguments.find("-.") != std::string::npos && lt < 0)
					{
						arguments.erase(arguments.find("-."),2);
						lt = 4;
					}
					else if(arguments.find("-") != std::string::npos && lt < 0)
					{
						arguments.erase(arguments.find("-"),1);
						lt = 1;
					}
					else if(arguments.find(":") != std::string::npos && lt < 0)
					{
						arguments.erase(arguments.find(":"),1);
						lt = 3;
					}
				
					if(arguments.find("+") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("+"),1);
						pt = 1;
					}
					else if(arguments.find("o") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("o"),1);
						pt = 6;
					}
					else if(arguments.find("*") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("*"),1);
						pt = 3;
					}
					else if(arguments.find(".") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("."),1);
						pt = 3;
					}
					else if(arguments.find("x") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("x"),1);
						pt = 2;
					}
					else if(arguments.find("*") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("*"),1);
						pt = 3;
					}
					else if(arguments.find("square") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("square"),6);
						pt = 4;
					}
					else if(arguments.find("s") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("s"),1);
						pt = 4;
					}
					else if(arguments.find("diamond") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("diamond"),7);
						pt = 12;
					}
					else if(arguments.find("d") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("d"),1);
						pt = 12;
					}
					else if(arguments.find("^") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("^"),1);
						pt = 8;
					}
					else if(arguments.find("v") != std::string::npos && pt < 0)
					{
						arguments.erase(arguments.find("v"),1);
						pt = 10;
					}
				
					if(arguments.find("k") != std::string::npos && color_str.empty())
					{
						arguments.erase(arguments.find("k"),1);
						color_str = "black";
					}
					else if(arguments.find("b") != std::string::npos && color_str.empty())
					{
						arguments.erase(arguments.find("b"),1);
						color_str = "blue";
					}
					else if(arguments.find("r") != std::string::npos && color_str.empty())
					{
						arguments.erase(arguments.find("r"),1);
						color_str = "red";
					}
					else if(arguments.find("g") != std::string::npos && color_str.empty())
					{
						arguments.erase(arguments.find("g"),1);
						color_str = "green";
					}
					else if(arguments.find("c") != std::string::npos && color_str.empty())
					{
						arguments.erase(arguments.find("c"),1);
						color_str = "cyan";
					}
					else if(arguments.find("m") != std::string::npos && color_str.empty())
					{
						arguments.erase(arguments.find("m"),1);
						color_str = "magenta";
					}
					else if(arguments.find("y") != std::string::npos && color_str.empty())
					{
						arguments.erase(arguments.find("y"),1);
						color_str = "yellow";
					}
				
					if(!arguments.empty())
					{
						throw FigureException("Unknown format string in Figure!");
					}
				}
			}
			else
			{
				lt = 1;
			}
			if(lt < 0 && pt < 0 && !color_str.empty())
			{
				lt = 1;
			}
			
			std::stringstream term_stream;
			term_stream << "wxt size 500,375 enhanced font 'Verdana,";
			term_stream << fontsize;
			term_stream << "' persist";
			g.set_terminal_std(term_stream.str());
			g.showonscreen();
			if(!hold_on_bool)
			{
				g.reset_plot();
			}
			g.cmd("set termoption dashed");
			g.cmd("set border linewidth 1.5");
			
			std::stringstream stream1;
			std::stringstream stream2;
			if(lt > 0 && pt > 0)
			{
				stream1 << "set style line " << num_plots << " linecolor rgb '";
				if(color_str.empty())
				{
					stream1 << "#";
					stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((num_plots-1) % colors.size(1),0));
					stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((num_plots-1) % colors.size(1),1));
					stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((num_plots-1) % colors.size(1),2));
				}
				else
				{
					stream1 << color_str;
				}
				stream1 << "' linetype " << lt << " lw " << lw;
				stream1 << " pt " << pt << " ps " << ps;
				g.cmd(stream1.str());
				
				stream2 << "linespoints ls " << num_plots;
			}
			else if(lt > 0)
			{
				stream1 << "set style line " << num_plots << " linecolor rgb '";
				if(color_str.empty())
				{
					stream1 << "#";
					stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((num_plots-1) % colors.size(1),0));
					stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((num_plots-1) % colors.size(1),1));
					stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((num_plots-1) % colors.size(1),2));
				}
				else
				{
					stream1 << color_str;
				}
				stream1 << "' linetype " << lt << " lw " << lw;
				g.cmd(stream1.str());
				
				stream2 << "lines ls " << num_plots;
			}
			else if(pt > 0)
			{
				stream2 << "points pt " << pt << " ps " << ps << " lc rgb '";
				if(color_str.empty())
				{
					stream2 << "#";
					stream2 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((num_plots-1) % colors.size(1),0));
					stream2 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((num_plots-1) % colors.size(1),1));
					stream2 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((num_plots-1) % colors.size(1),2));
				}
				else
				{
					stream2 << color_str;
				}
				stream2 << "'";
			}
			else
			{
				throw FigureException("Unknown error while plotting in Figure!");
			}
			g.set_style(stream2.str());
			g.plot_xy(x,y,legend_entry);
		}
		catch(GnuplotException ge)
		{
			std::cout << ge.what() << std::endl;
		}
		
		if(num_plots > x_plot_data.size())
		{
			x_plot_data.push_back(x);
			y_plot_data.push_back(y);
			plot_format.push_back(format);
			plot_linewidth.push_back(lw);
			plot_markersize.push_back(ps);
			plot_val.push_back(val);
		}
		
		return;
	}
	
	inline void Figure::xlabel(std::string xlabel_text)
	{
		try
		{
			g.set_xlabel(xlabel_text);
			g.replot();
		}
		catch(GnuplotException ge)
		{
			std::cout << ge.what() << std::endl;
		}
		m_xlabel = xlabel_text;
		return;
	}
	
	inline void Figure::ylabel(std::string ylabel_text)
	{
		try
		{
			g.set_ylabel(ylabel_text);
			g.replot();
		}
		catch(GnuplotException ge)
		{
			std::cout << ge.what() << std::endl;
		}
		m_ylabel = ylabel_text;
		return;
	}
	
	inline void Figure::title(std::string title_text)
	{
		try
		{
			g.set_title(title_text);
			g.replot();
		}
		catch(GnuplotException ge)
		{
			std::cout << ge.what() << std::endl;
		}
		m_title = title_text;
		return;
	}
	
	inline void Figure::grid_on()
	{
		try
		{
			g.set_grid();
			g.replot();
		}
		catch(GnuplotException ge)
		{
			std::cout << ge.what() << std::endl;
		}
		grid_on_bool = true;
		return;
	}
	
	inline void Figure::grid_off()
	{
		try
		{
			g.unset_grid();
			g.replot();
		}
		catch(GnuplotException ge)
		{
			std::cout << ge.what() << std::endl;
		}
		grid_on_bool = false;
		return;
	}
	
	inline void Figure::hold_on()
	{
		hold_on_bool = true;
		return;
	}
	
	inline void Figure::hold_off()
	{
		hold_on_bool = false;
		return;
	}
	
	inline void Figure::legend(const std::initializer_list<std::string> lst)
	{
		if(lst.size() <= 0)
		{
			throw FigureException("Cannot create empty legend!");
		}
		if(lst.size() > num_plots)
		{
			throw FigureException("Error! You tried to create a legend with more entries than plots!");
		}
		int ii = 0;
		legend_entries = std::vector<std::string>(lst.size());
		for(const auto& l : lst)
		{
			legend_entries[ii] = l;
			ii++;
		}
		replot_all();
	}
	
	
	inline void Figure::replot_all()
	{
		num_plots = 0;
		try
		{
			g.reset_all();
		}
		catch(GnuplotException ge)
		{
			std::cout << ge.what() << std::endl;
		}

		for(int jj = 0; jj < x_plot_data.size(); jj++)
		{
			if(legend_entries.size() > jj)
			{
				plot(x_plot_data[jj],y_plot_data[jj],plot_format[jj],plot_val[jj],plot_linewidth[jj],plot_markersize[jj],legend_entries[jj]);
			}
			else
			{
				plot(x_plot_data[jj],y_plot_data[jj],plot_format[jj],plot_val[jj],plot_linewidth[jj],plot_markersize[jj],"");
			}
		}
		
		if(!m_xlabel.empty())
		{
			xlabel(m_xlabel);
		}
		if(!m_ylabel.empty())
		{
			ylabel(m_ylabel);
		}
		if(!m_title.empty())
		{
			title(m_title);
		}
		if(grid_on_bool)
		{
			grid_on();
		}
	}
	
	inline void Figure::set(std::string property, double val)
	{
		std::transform(property.begin(), property.end(), property.begin(), ::tolower);
		if(property.compare("fontsize") == 0)
		{
			fontsize = (int)round(val);
			replot_all();
		}
		else
		{
			throw FigureException("Unknown property in set!");
		}
	}

}

#endif
