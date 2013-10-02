#ifndef vector_k_KEYCPP_H_
#define vector_k_KEYCPP_H_

#include <vector>

namespace keycpp
{
    template <class T>
    class  vector_k
    {
    public:

        typedef T *iterator;
        typedef const T* const_iterator;

        vector_k();
        vector_k(unsigned int size);
        vector_k(unsigned int size, const T &initial);
        vector_k(T *ptr, unsigned int size, unsigned int pinc);
		vector_k(const std::initializer_list<T>& lst);
        vector_k(const vector_k<T> &v);      
        vector_k(std::vector<T> const &v);      
        ~vector_k();

        unsigned int capacity() const;
        unsigned int size() const;
        bool empty() const;
        iterator begin();
        iterator end();
        const_iterator begin() const;
        const_iterator end() const;
        T& front();
        T& back();
        void push_back(const T &value); 
        void pop_back();  

        void reserve(unsigned int capacity);   
        void resize(unsigned int size);   

        T& operator[](const unsigned int index);  
        const T& operator[](const unsigned int index) const;
        vector_k<T>& operator=(const vector_k<T>& v);
        vector_k<T>& operator=(const std::vector<T>& v);
        void clear();
        
        operator std::vector<T>()
        {
            std::vector<T> v1(my_size);
            for(unsigned int ii = 0; ii < my_size; ii++)
            {
                v1[ii] = buffer[ii*inc];
            }
            return v1;
        };

    private:
        unsigned int my_size;
        unsigned int my_capacity;
        unsigned int inc = 1;
        T *buffer;
        bool dontFree = true;
    };

    // Your code goes here ...
    template<class T>
    vector_k<T>::vector_k()
    {
        my_capacity = 0;
        my_size = 0;
        buffer = 0;
        dontFree = true;
    }

    template<class T>
    vector_k<T>::vector_k(const vector_k<T> &v)
    {
        my_size = v.my_size;
        my_capacity = v.my_capacity;
        dontFree = false;
        inc = 1;
        buffer = new T[my_size];  
        for(unsigned int ii = 0; ii < my_size; ii++)
        {
            buffer[ii*inc] = v.buffer[ii];  
        }
    }

    template<class T>
    vector_k<T>::vector_k(unsigned int size)
    {
        my_capacity = size;
        my_size = size;
        dontFree = false;
        inc = 1;
        buffer = new T[size];
    }

    template<class T>
    vector_k<T>::vector_k(unsigned int size, const T &initial)
    {
        my_size = size;
        my_capacity = size;
        dontFree = false;
        inc = 1;
        buffer = new T[size];
        for(unsigned int ii = 0; ii < size; ii++)
        {
            buffer[ii*inc] = initial;
        }
    }

    template<class T>
    vector_k<T>::vector_k(T *ptr, unsigned int size, unsigned int pinc)
    {
        my_capacity = size;
        my_size = size;
        inc = pinc;
        dontFree = true;
        buffer = ptr;
    }
    
	template<class T>
	vector_k<T>::vector_k(const std::initializer_list<T>& lst) : vector_k(lst.size())
	{
		int ii = 0;
		for(const auto& l : lst)
		{
			buffer[ii*inc] = l;
			ii++;
		}
	}

    template<class T>
    vector_k<T>& vector_k<T>::operator=(const vector_k<T> &v)
    {
        if(my_size != v.size())
        {
            resize(v.size());
        }
        for(unsigned int ii = 0; ii < my_size; ii++)
        {
            buffer[ii*inc] = v.buffer[ii];
        }
        
        return *this;
    }

    template<class T>
    vector_k<T>& vector_k<T>::operator=(std::vector<T> const &v)
    {
        if(my_size != v.size())
        {
            resize(v.size());
        }
        for(unsigned int ii = 0; ii < my_size; ii++)
        {
            buffer[ii*inc] = v[ii];
        }
        
        return *this;
    }

    template<class T>
    typename vector_k<T>::iterator vector_k<T>::begin()
    {
        return buffer;
    }

    template<class T>
    typename vector_k<T>::iterator vector_k<T>::end()
    {
        return buffer + size();
    }

    template<class T>
    typename vector_k<T>::const_iterator vector_k<T>::begin() const
    {
        return buffer;
    }

    template<class T>
    typename vector_k<T>::const_iterator vector_k<T>::end() const
    {
        return buffer + size();
    }

    template<class T>
    T& vector_k<T>::front()
    {
        return buffer[0*inc];
    }

    template<class T>
    T& vector_k<T>::back()
    {
        return buffer[(my_size-1)*inc];
    }

    template<class T>
    void vector_k<T>::push_back(const T& v)
    {
        if(my_size >= my_capacity)
        {
            reserve(2*my_size + 5);
        }
        buffer[(my_size++)*inc] = v;
    }

    template<class T>
    void vector_k<T>::pop_back()
    {
        my_size--;
    }

    template<class T>
    void vector_k<T>::reserve(unsigned int capacity)
    {
        if(buffer == 0)
        {
            my_size = 0;
            my_capacity = 0;
        }
        T *Newbuffer = new T[capacity];
        unsigned int l_Size = (capacity < my_size)? capacity : my_size;

        for(unsigned int ii = 0; ii < l_Size; ii++)
        {
            Newbuffer[ii] = buffer[ii*inc];
        }
        inc = 1;

        my_capacity = capacity;
        if(!dontFree)
        {
            delete[] buffer;
            dontFree = true;
        }
        dontFree = false;
        buffer = Newbuffer;
    }

    template<class T>
    unsigned int vector_k<T>::size() const
    {
        return my_size;
    }

    template<class T>
    bool vector_k<T>::empty() const
    {
        bool val = false;
        if(my_size == 0)
        {
            val = true;
        }
        return val;
    }

    template<class T>
    void vector_k<T>::resize(unsigned int size)
    {
        reserve(size);
        my_size = size;
    }

    template<class T>
    T& vector_k<T>::operator[](const unsigned int index)
    {
        return buffer[index*inc];
    }  
    
    template<class T>
    const T& vector_k<T>::operator[](const unsigned int index) const
    {
        return buffer[index*inc];
    }  

    template<class T>
    unsigned int vector_k<T>::capacity() const
    {
        return my_capacity;
    }

    template<class T>
    vector_k<T>::~vector_k()
    {
        if(!dontFree)
        {
            delete[] buffer;
        }
    }
    template <class T>
    void vector_k<T>::clear()
    {
        my_capacity = 0;
        my_size = 0;
        if(!dontFree)
        {
            delete[] buffer;
        }
        buffer = 0;
    }
}

#endif
