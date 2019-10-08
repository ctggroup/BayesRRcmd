branch=$(git rev-parse --abbrev-ref HEAD)
ver=$(git describe --always) 
echo "
#ifndef VERSION_H
#define VERSION_H
extern char const *const GIT_COMMIT=\"Branch:$branch, commit:$ver\";
#endif
" > src/version.h 
