from distutils.core import setup, Extension
import sysconfig

def main():
    CFLAGS = ['-g', '-Wall', '-std=c99', '-fopenmp', '-mavx', '-mfma', '-pthread', '-O3']
    LDFLAGS = ['-fopenmp']
    # Use the setup function we imported and set up the modules.
    # You may find this reference helpful: https://docs.python.org/3.6/extending/building.html
    
    module = Extension('numc',
                       extra_compile_args=CFLAGS,
                       extra_link_args=LDFLAGS,
                       sources=['src/matrix.c','src/numc.c'])

    setup(name='numc',
          ext_modules=[module])

    

if __name__ == "__main__":
    main()
