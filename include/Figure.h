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


/* TODO:
add function to resize or set the size of the figure
add semilogx, semilogy, loglog
add subplot capability
add plotyy (low priority)
add contour, contourf, ezcontour, ezcontourf
add ezplot
add quiver
add scatter
add plot3
add surface
add scatter3
*/

namespace keycpp
{
	class FigureException : public std::runtime_error
	{
		public:
			FigureException(const std::string &msg) : std::runtime_error(msg){}
	};
	
	class Plots
	{
	public:
	    Plots();
		bool hold_on_bool = false;
		int num_plots = 0;
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
		double ymin;
		double ymax;
		double xmin;
		double xmax;
		bool grid_on_bool = false;
		bool logscale_x = false;
		bool logscale_y = false;
        std::ostringstream cmdstr;
	};
	
	inline Plots::Plots()
	{
		ymin = nan("");
		ymax = nan("");
		xmin = nan("");
		xmax = nan("");
	}
	
	class Figure
	{
	private:
		Gnuplot g;
		matrix<double> colors;
		std::vector<Plots> p;
		int fontsize = 10;
		bool multiplot = false;
		int current_plot = 0;
		int multi_rows = -1;
		int multi_cols = -1;
		bool final_replot = false;
		int m_width = 500;
		int m_height = 375;
	
	public:
		Figure();
		~Figure();
		template<class U, class T> void plot(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1);
		template<class U, class T> void plot(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2);
		template<class U, class T> void plot(std::vector<U> x, std::vector<T> y, std::string arguments = "", double val = -1, double lw = 2, double ps = 1.5, std::string legend_entry = "");
		template<class T> void plot(std::vector<T> y, std::string format, std::string property1, double val1);
		template<class T> void plot(std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2);
		template<class T> void plot(std::vector<T> y, std::string arguments = "", double val = -1, double lw = 2, double ps = 1.5, std::string legend_entry = "");
		template<class U, class T> void semilogx(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1);
		template<class U, class T> void semilogx(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2);
		template<class U, class T> void semilogx(std::vector<U> x, std::vector<T> y, std::string arguments = "", double val = -1, double lw = 2, double ps = 1.5, std::string legend_entry = "");
		template<class U, class T> void semilogy(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1);
		template<class U, class T> void semilogy(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2);
		template<class U, class T> void semilogy(std::vector<U> x, std::vector<T> y, std::string arguments = "", double val = -1, double lw = 2, double ps = 1.5, std::string legend_entry = "");
		template<class U, class T> void loglog(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1);
		template<class U, class T> void loglog(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2);
		template<class U, class T> void loglog(std::vector<U> x, std::vector<T> y, std::string arguments = "", double val = -1, double lw = 2, double ps = 1.5, std::string legend_entry = "");
		void xlabel(std::string xlabel_text);
		void ylabel(std::string ylabel_text);
		void title(std::string title_text);
		void grid_on();
		void grid_off();
		void hold_on();
		void hold_off();
		void subplot(int mrows, int mcols, int index);
		void legend(std::initializer_list<std::string> lst);
		void ylim(std::initializer_list<double> lst);
		void xlim(std::initializer_list<double> lst);
		void replot_all();
		void setFontsize(int p_fontsize) {fontsize = p_fontsize;};
		int getFontsize() {return fontsize;};
		void set(std::string property, double val);
		void set(std::string property, std::initializer_list<int> list);
	};
	
	inline Figure::Figure() try : g("lines")
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
		
		p = std::vector<Plots>(1);
	}
    catch(GnuplotException ge)
    {
	    std::cout << ge.what() << std::endl;
    }
	
	inline Figure::~Figure() 
	{
	    final_replot = true;
	    replot_all();
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
		std::transform(arguments.begin(), arguments.end(), arguments.begin(), ::tolower);
		format = arguments;
		int lt = -1;
		int pt = -1;
		std::string color_str = "";
		
		p[current_plot].num_plots++;
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
		
		if(final_replot)
		{
		    try
		    {
			    if(current_plot == 0 && p[current_plot].num_plots == 1)
			    {
			        std::stringstream term_stream;
			        term_stream << "wxt size ";
			        term_stream << m_width << "," << m_height << " enhanced font 'Verdana,";
			        term_stream << fontsize;
			        term_stream << "' persist";
			        g.set_terminal_std(term_stream.str());
			        g.showonscreen();
			        if(multiplot)
			        {
			            std::stringstream multi_stream;
			            multi_stream << "set multiplot layout ";
			            multi_stream << multi_rows << "," << multi_cols;
			            multi_stream << " columnsfirst";
			            g.cmd(multi_stream.str());
			        }
			        g.cmd("set termoption dashed");
			        g.cmd("set border linewidth 1.5");
			    }
			
			    if(!p[current_plot].hold_on_bool)
			    {
				    g.reset_plot();
			    }
			
			    if(p[current_plot].xmin == p[current_plot].xmin && p[current_plot].xmax == p[current_plot].xmax) // Check for NaNs
			    {
				    std::stringstream temp_stream;
				    temp_stream << "set xrange [";
				    temp_stream << p[current_plot].xmin << ":" << p[current_plot].xmax << "]";
				    g.cmd(temp_stream.str());
			    }
			    if(p[current_plot].ymin == p[current_plot].ymin && p[current_plot].ymax == p[current_plot].ymax) // Check for NaNs
			    {
				    std::stringstream temp_stream;
				    temp_stream << "set yrange [";
				    temp_stream << p[current_plot].ymin << ":" << p[current_plot].ymax << "]";
				    g.cmd(temp_stream.str());
			    }
			
			    std::stringstream stream1;
			    std::stringstream stream2;
			    if(lt > 0 && pt > 0)
			    {
				    stream1 << "set style line " << p[current_plot].num_plots << " linecolor rgb '";
				    if(color_str.empty())
				    {
					    stream1 << "#";
					    stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((p[current_plot].num_plots-1) % colors.size(1),0));
					    stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((p[current_plot].num_plots-1) % colors.size(1),1));
					    stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((p[current_plot].num_plots-1) % colors.size(1),2));
				    }
				    else
				    {
					    stream1 << color_str;
				    }
				    stream1 << "' linetype " << lt << " lw " << lw;
				    stream1 << " pt " << pt << " ps " << ps;
				    g.cmd(stream1.str());
				
				    stream2 << "linespoints ls " << p[current_plot].num_plots;
			    }
			    else if(lt > 0)
			    {
				    stream1 << "set style line " << p[current_plot].num_plots << " linecolor rgb '";
				    if(color_str.empty())
				    {
					    stream1 << "#";
					    stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((p[current_plot].num_plots-1) % colors.size(1),0));
					    stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((p[current_plot].num_plots-1) % colors.size(1),1));
					    stream1 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((p[current_plot].num_plots-1) % colors.size(1),2));
				    }
				    else
				    {
					    stream1 << color_str;
				    }
				    stream1 << "' linetype " << lt << " lw " << lw;
				    g.cmd(stream1.str());
				
				    stream2 << "lines ls " << p[current_plot].num_plots;
			    }
			    else if(pt > 0)
			    {
				    stream2 << "points pt " << pt << " ps " << ps << " lc rgb '";
				    if(color_str.empty())
				    {
					    stream2 << "#";
					    stream2 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((p[current_plot].num_plots-1) % colors.size(1),0));
					    stream2 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((p[current_plot].num_plots-1) % colors.size(1),1));
					    stream2 << std::setfill ('0') << std::setw(2) << std::hex << (int)round(255*colors((p[current_plot].num_plots-1) % colors.size(1),2));
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
			
			
			    if(!p[current_plot].m_xlabel.empty())
			    {
			        g.set_xlabel(p[current_plot].m_xlabel);
			    }
			    if(!p[current_plot].m_ylabel.empty())
			    {
			        g.set_ylabel(p[current_plot].m_ylabel);
			    }
			    if(!p[current_plot].m_title.empty())
			    {
			        g.set_title(p[current_plot].m_title);
			    }
			    if(p[current_plot].grid_on_bool)
			    {
			        g.set_grid();
			    }
			    if(p[current_plot].logscale_x)
			    {
			        g.set_xlogscale();
			    }
			    if(p[current_plot].logscale_y)
			    {
			        g.set_ylogscale();
			    }
			
			    if(multiplot)
			    {
                    std::ofstream tmp;
                    std::string name = g.create_tmpfile(tmp);
                    if(name == "")
                    {
                        throw FigureException("Error creating temporary file!");
                    }

                    for(int ii = 0; ii < x.size(); ii++)
                    {
                        tmp << x[ii] << " " << y[ii] << std::endl;
                    }
                    tmp.flush();
                    tmp.close();
                    
                    if(p[current_plot].num_plots > 1)
                    {
                        p[current_plot].cmdstr << ", ";
                    }
                    else
                    {
                        p[current_plot].cmdstr << "plot ";
                    }

                    p[current_plot].cmdstr << "\"" << name << "\" using 1:2";

                    if(legend_entry.empty())
                    {
                        p[current_plot].cmdstr << " notitle ";
                    }
                    else
                    {
                        p[current_plot].cmdstr << " title \"" << legend_entry << "\" ";
                    }
                    p[current_plot].cmdstr << "with " << g.get_style();
			    }
			    else
			    {
			        g.plot_xy(x,y,legend_entry);
			    }
			
			
		        if(multiplot && current_plot == p.size()-1 && p[current_plot].num_plots == p[current_plot].x_plot_data.size())
		        {
	                if(p[current_plot].cmdstr.str().length() != 0)
	                {
                        g.cmd(p[current_plot].cmdstr.str());
	                }
		            g.cmd("unset multiplot");
		        }
		    }
		    catch(GnuplotException ge)
		    {
			    std::cout << ge.what() << std::endl;
		    }
		}
		
		if(p[current_plot].num_plots > p[current_plot].x_plot_data.size())
		{
			p[current_plot].x_plot_data.push_back(x);
			p[current_plot].y_plot_data.push_back(y);
			p[current_plot].plot_format.push_back(format);
			p[current_plot].plot_linewidth.push_back(lw);
			p[current_plot].plot_markersize.push_back(ps);
			p[current_plot].plot_val.push_back(val);
		}
		
		return;
	}
	
	template<class T> void Figure::plot(std::vector<T> y, std::string format, std::string property1, double val1)
	{
	    std::vector<T> x(y.size());
	    for(int ii = 0; ii < y.size(); ii++)
	    {
	        x[ii] = ii;
	    }
	    plot(x,y,format,property1,val1);
	}
	
	template<class T> void Figure::plot(std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2)
	{
	    std::vector<T> x(y.size());
	    for(int ii = 0; ii < y.size(); ii++)
	    {
	        x[ii] = ii;
	    }
	    plot(x,y,format,property1,val1,property2,val2);
	}
	
	template<class T> void Figure::plot(std::vector<T> y, std::string arguments, double val, double lw, double ps, std::string legend_entry)
	{
	    std::vector<T> x(y.size());
	    for(int ii = 0; ii < y.size(); ii++)
	    {
	        x[ii] = ii;
	    }
	    plot(x,y,arguments,val,lw,ps,legend_entry);
	}
	
	template<class U, class T> void Figure::semilogx(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1)
	{
	    p[current_plot].logscale_x = true;
	    plot(x,y,format,property1,val1);
	}
	
	template<class U, class T> void Figure::semilogx(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2)
	{
	    p[current_plot].logscale_x = true;
	    plot(x,y,format,property1,val1,property2,val2);
	}
	
	template<class U, class T> void Figure::semilogx(std::vector<U> x, std::vector<T> y, std::string arguments, double val, double lw, double ps, std::string legend_entry)
	{
	    p[current_plot].logscale_x = true;
	    plot(x,y,arguments,val,lw,ps,legend_entry);
	}
	
	template<class U, class T> void Figure::semilogy(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1)
	{
	    p[current_plot].logscale_y = true;
	    plot(x,y,format,property1,val1);
	}
	
	template<class U, class T> void Figure::semilogy(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2)
	{
	    p[current_plot].logscale_y = true;
	    plot(x,y,format,property1,val1,property2,val2);
	}
	
	template<class U, class T> void Figure::semilogy(std::vector<U> x, std::vector<T> y, std::string arguments, double val, double lw, double ps, std::string legend_entry)
	{
	    p[current_plot].logscale_y = true;
	    plot(x,y,arguments,val,lw,ps,legend_entry);
	}
	
	template<class U, class T> void Figure::loglog(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1)
	{
	    p[current_plot].logscale_x = true;
	    p[current_plot].logscale_y = true;
	    plot(x,y,format,property1,val1);
	}
	
	template<class U, class T> void Figure::loglog(std::vector<U> x, std::vector<T> y, std::string format, std::string property1, double val1, std::string property2, double val2)
	{
	    p[current_plot].logscale_x = true;
	    p[current_plot].logscale_y = true;
	    plot(x,y,format,property1,val1,property2,val2);
	}
	
	template<class U, class T> void Figure::loglog(std::vector<U> x, std::vector<T> y, std::string arguments, double val, double lw, double ps, std::string legend_entry)
	{
	    p[current_plot].logscale_x = true;
	    p[current_plot].logscale_y = true;
	    plot(x,y,arguments,val,lw,ps,legend_entry);
	}
	
	inline void Figure::xlabel(std::string xlabel_text)
	{
		p[current_plot].m_xlabel = xlabel_text;
		return;
	}
	
	inline void Figure::ylabel(std::string ylabel_text)
	{
		p[current_plot].m_ylabel = ylabel_text;
		return;
	}
	
	inline void Figure::title(std::string title_text)
	{
		p[current_plot].m_title = title_text;
		return;
	}
	
	inline void Figure::grid_on()
	{
		p[current_plot].grid_on_bool = true;
		return;
	}
	
	inline void Figure::grid_off()
	{
		p[current_plot].grid_on_bool = false;
		return;
	}
	
	inline void Figure::hold_on()
	{
		p[current_plot].hold_on_bool = true;
		return;
	}
	
	inline void Figure::hold_off()
	{
		p[current_plot].hold_on_bool = false;
		return;
	}
	
	inline void Figure::legend(const std::initializer_list<std::string> lst)
	{
		if(lst.size() <= 0)
		{
			throw FigureException("Cannot create empty legend!");
		}
		if(lst.size() > p[current_plot].num_plots)
		{
			throw FigureException("Error! You tried to create a legend with more entries than plots!");
		}
		int ii = 0;
		p[current_plot].legend_entries = std::vector<std::string>(lst.size());
		for(const auto& l : lst)
		{
			p[current_plot].legend_entries[ii] = l;
			ii++;
		}
	}
	
	
	inline void Figure::replot_all()
	{
	    int N = current_plot;
	    for(current_plot = 0; current_plot <= N; current_plot++)
	    {
		    p[current_plot].num_plots = 0;

		    for(int jj = 0; jj < p[current_plot].x_plot_data.size(); jj++)
		    {
			    if(p[current_plot].legend_entries.size() > jj)
			    {
				    plot(p[current_plot].x_plot_data[jj],p[current_plot].y_plot_data[jj],p[current_plot].plot_format[jj],p[current_plot].plot_val[jj],p[current_plot].plot_linewidth[jj],p[current_plot].plot_markersize[jj],p[current_plot].legend_entries[jj]);
			    }
			    else
			    {
				    plot(p[current_plot].x_plot_data[jj],p[current_plot].y_plot_data[jj],p[current_plot].plot_format[jj],p[current_plot].plot_val[jj],p[current_plot].plot_linewidth[jj],p[current_plot].plot_markersize[jj],"");
			    }
		    }
		    
		    if(!p[current_plot].cmdstr.str().empty() && current_plot != p.size()-1)
		    {
                g.cmd(p[current_plot].cmdstr.str());
		    }
	    }
	    current_plot = N;
	}
	
	inline void Figure::set(std::string property, double val)
	{
		std::transform(property.begin(), property.end(), property.begin(), ::tolower);
		if(property.compare("fontsize") == 0)
		{
			fontsize = (int)round(val);
		}
		else
		{
			throw FigureException("Unknown property in set!");
		}
	}
	
	inline void Figure::set(std::string property, std::initializer_list<int> list)
	{
		std::transform(property.begin(), property.end(), property.begin(), ::tolower);
		if(property.compare("position") == 0)
		{
		    // Note for compatibility with MATLAB, I have kept the position values
		    // but we only use the height and width values.
		    int ii = 0;
		    for(const auto& l : list)
		    {
			    if(ii == 2)
	            {
			        m_width = l;
		        }
		        else if(ii == 3)
	            {
		            m_height = l;
		        }
		        else if(ii > 3)
		        {
		            throw FigureException("Too many values while trying to set position!");
		        }
			    ii++;
		    }
		}
		else
		{
			throw FigureException("Unknown property in set!");
		}
	}
	
	inline void Figure::ylim(const std::initializer_list<double> lst)
	{
		if(lst.size() <= 0)
		{
			throw FigureException("ylim() called with no limits!");
		}
		if(lst.size() > 2)
		{
			throw FigureException("Error! More than 2 limits provided!");
		}
		int ii = 0;
		for(const auto& l : lst)
		{
			if(ii == 0)
			{
				p[current_plot].ymin = l;
			}
			else if(ii == 1)
			{
				p[current_plot].ymax = l;
			}
			ii++;
		}
	}
	
	inline void Figure::xlim(const std::initializer_list<double> lst)
	{
		if(lst.size() <= 0)
		{
			throw FigureException("ylim() called with no limits!");
		}
		if(lst.size() > 2)
		{
			throw FigureException("Error! More than 2 limits provided!");
		}
		int ii = 0;
		for(const auto& l : lst)
		{
			if(ii == 0)
			{
				p[current_plot].xmin = l;
			}
			else if(ii == 1)
			{
				p[current_plot].xmax = l;
			}
			ii++;
		}
	}
	
	inline void Figure::subplot(int mrows, int mcols, int index)
	{
	    if((mrows != multi_rows && multi_rows > 0) || (mcols != multi_cols && multi_cols > 0))
	    {
	        throw FigureException("Layout of subplot does not match previous calls!");
	    }
	    if(mrows <= 0 || mcols <= 0 || index < 0 || index >= mrows*mcols)
	    {
	        throw FigureException("Illegal argument in subplot!");
	    }
	    multiplot = true;
	    multi_rows = mrows;
	    multi_cols = mcols;
	    current_plot = index;
	    if(p.size() != mrows*mcols)
	    {
	        p = std::vector<Plots>(mrows*mcols);
	    }
	}
}

#endif
