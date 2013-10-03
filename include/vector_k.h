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
        typedef size_t size_type;

        vector_k();
        vector_k(size_t size);
        vector_k(size_t size, const T &initial);
        vector_k(T *ptr, size_t size, size_t pinc);
		vector_k(const std::initializer_list<T>& lst);
        vector_k(const vector_k<T> &v);      
        vector_k(const std::vector<T> &v);      
        ~vector_k();

        size_t capacity() const;
        size_t size() const;
        bool empty() const;
        iterator begin();
        iterator end();
        const_iterator begin() const;
        const_iterator end() const;
        T& front();
        T& back();
        void push_back(const T &value); 
        void pop_back();  

        void reserve(size_t capacity);   
        void resize(size_t size);   

        T& operator[](const size_t index);  
        const T& operator[](const size_t index) const;
        vector_k<T>& operator=(const vector_k<T>& v);
        vector_k<T>& operator=(const std::vector<T>& v);
        void clear();
        
        operator std::vector<T>()
        {
            std::vector<T> v1(my_size);
            for(size_t ii = 0; ii < my_size; ii++)
            {
                v1[ii] = buffer[ii*inc];
            }
            return v1;
        };

    private:
        size_t my_size;
        size_t my_capacity;
        size_t inc = 1;
        T *buffer;
        bool dontFree = true;
    };

    // Your code goes here ...
    template<class T>
    vector_k<T>::vector_k() : my_size(0), my_capacity(0), buffer(nullptr)
    {}

    template<class T>
    vector_k<T>::vector_k(const vector_k<T> &v) : my_size(v.my_size), my_capacity(v.my_capacity), buffer(new T[v.my_size])
    {
        dontFree = false;
        inc = 1;
        for(size_t ii = 0; ii < my_size; ii++)
        {
            buffer[ii*inc] = v.buffer[ii];  
        }
    }

    template<class T>
    vector_k<T>::vector_k(const std::vector<T> &v) : my_size(v.size()), my_capacity(v.capacity()), buffer(new T[v.size()])
    {
        dontFree = false;
        inc = 1;
        for(size_t ii = 0; ii < my_size; ii++)
        {
            buffer[ii*inc] = v[ii];  
        }
    }

    template<class T>
    vector_k<T>::vector_k(size_t p_size) : my_size(p_size), my_capacity(p_size), buffer(new T[p_size])
    {
        dontFree = false;
        inc = 1;
    }

    template<class T>
    vector_k<T>::vector_k(size_t p_size, const T &initial) : my_size(p_size), my_capacity(p_size), buffer(new T[p_size])
    {
        dontFree = false;
        inc = 1;
        for(size_t ii = 0; ii < p_size; ii++)
        {
            buffer[ii*inc] = initial;
        }
    }

    template<class T>
    vector_k<T>::vector_k(T *ptr, size_t p_size, size_t pinc) : my_size(p_size), my_capacity(p_size), buffer(nullptr)
    {
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
        for(size_t ii = 0; ii < my_size; ii++)
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
        for(size_t ii = 0; ii < my_size; ii++)
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
    void vector_k<T>::reserve(size_t p_capacity)
    {
        if(buffer == nullptr)
        {
            my_size = 0;
            my_capacity = 0;
        }
        T *Newbuffer = new T[p_capacity];
        size_t l_Size = (p_capacity < my_size)? p_capacity : my_size;

        for(size_t ii = 0; ii < l_Size; ii++)
        {
            Newbuffer[ii] = buffer[ii*inc];
        }
        inc = 1;

        my_capacity = p_capacity;
        if(!dontFree)
        {
            delete[] buffer;
            dontFree = true;
        }
        dontFree = false;
        buffer = Newbuffer;
    }

    template<class T>
    size_t vector_k<T>::size() const
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
    void vector_k<T>::resize(size_t p_size)
    {
        reserve(p_size);
        my_size = p_size;
    }

    template<class T>
    T& vector_k<T>::operator[](const size_t index)
    {
        return buffer[index*inc];
    }  
    
    template<class T>
    const T& vector_k<T>::operator[](const size_t index) const
    {
        return buffer[index*inc];
    }  

    template<class T>
    size_t vector_k<T>::capacity() const
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
        buffer = nullptr;
    }
}

#endif
