import os
import multiprocessing

SetOption('num_jobs', multiprocessing.cpu_count())

pkgs = [ ]

env = Environment()
env['ENV']['TERM'] = os.environ['TERM']
env['CXXCOMSTR'] = " (CXX) $TARGET"
env['ARCOMSTR'] = " (AR)  $TARGET"
env['RANLIBCOMSTR'] = " (RAN) $TARGET"
env['LINKCOMSTR'] = " (LD)  $TARGET"
if len(pkgs) != 0:
    env.ParseConfig('pkg-config --cflags --libs '+" ".join(pkgs))
env.ParseConfig('getfem-config --cflags --libs')

env.Append(CXXFLAGS = ['-pedantic', '-Wall', '-Wextra', '-Wno-unused-parameter' ]) #, '-fmax-errors=5'])
env.Append(LINKFLAGS = [ ])

if int(ARGUMENTS.get('release', 0)) == 1:
    env.Append(CXXFLAGS = ['-O3'])
else:
    env.Append(CXXFLAGS = ['-ggdb', '-O', '-rdynamic'])
    env.Append(LINKFLAGS = ['-rdynamic'])

env.Append(CPPPATH = ['#src','#../common'])

env.VariantDir('build', 'src', duplicate=0)
env.SConscript(['build/SConscript' ], 'env')
