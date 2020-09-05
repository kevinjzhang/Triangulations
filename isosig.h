#ifndef ISO_SIG_H
#define ISO_SIG_H

#include <vector>
#include "information.h"

#include<triangulation/dim3.h>
#include<triangulation/example3.h>
#include<triangulation/detail/triangulation.h>
#include<triangulation/detail/isosig-impl.h>

using namespace regina;
using namespace detail;

class IsoSig {
private:
    IsoSig();

    //This is a copy of the function from regina
    template <int dim>
    static std::string isoSigFrom (Triangulation<3>* triangulation, size_t simp, 
            const Perm<dim+1>& vertices, Isomorphism<dim>* relabelling) {
        size_t nSimp = triangulation->size();
        size_t nFacets = ((dim + 1) * triangulation->size() + triangulation->countBoundaryFacets()) / 2;
        char* facetAction = new char[nFacets];
        size_t* joinDest = new size_t[nFacets];
        typedef typename Perm<dim+1>::Index PermIndex;
        PermIndex* joinGluing = new PermIndex[nFacets];
        ptrdiff_t* image = new ptrdiff_t[nSimp];
        Perm<dim+1>* vertexMap = new Perm<dim+1>[nSimp];

        // The preimage for each simplex:
        ptrdiff_t* preImage = new ptrdiff_t[nSimp];

        // ---------------------------------------------------------------------
        // Looping variables
        // ---------------------------------------------------------------------
        size_t facetPos, joinPos, nextUnusedSimp;
        size_t simpImg, simpSrc, dest;
        unsigned facetImg, facetSrc;
        const Simplex<dim>* s;

        // ---------------------------------------------------------------------
        // The code!
        // ---------------------------------------------------------------------

        std::fill(image, image + nSimp, -1);
        std::fill(preImage, preImage + nSimp, -1);

        image[simp] = 0;
        vertexMap[simp] = vertices.inverse();
        preImage[0] = simp;

        facetPos = 0;
        joinPos = 0;
        nextUnusedSimp = 1;

        // To obtain a canonical isomorphism, we must run through the simplices
        // and their facets in image order, not preimage order.
        //
        // This main loop is guaranteed to exit when (and only when) we have
        // exhausted a single connected component of the triangulation.
        for (simpImg = 0; simpImg < nSimp && preImage[simpImg] >= 0; ++simpImg) {
            simpSrc = preImage[simpImg];
            s = triangulation->simplex(simpSrc);

            for (facetImg = 0; facetImg <= dim; ++facetImg) {
                facetSrc = vertexMap[simpSrc].preImageOf(facetImg);

                // INVARIANTS (held while we stay within a single component):
                // - nextUnusedSimp > simpImg
                // - image[simpSrc], preImage[image[simpSrc]] and vertexMap[simpSrc]
                //   are already filled in.

                // Work out what happens to our source facet.
                if (! s->adjacentSimplex(facetSrc)) {
                    // A boundary facet.
                    facetAction[facetPos++] = 0;
                    continue;
                }

                // We have a real gluing.  Is it a gluing we've already seen
                // from the other side?
                dest = s->adjacentSimplex(facetSrc)->index();

                if (image[dest] >= 0)
                    if (image[dest] < image[simpSrc] ||
                            (dest == simpSrc &&
                            vertexMap[simpSrc][s->adjacentFacet(facetSrc)]
                            < vertexMap[simpSrc][facetSrc])) {
                        // Yes.  Just skip this gluing entirely.
                        continue;
                    }

                // Is it a completely new simplex?
                if (image[dest] < 0) {
                    // Yes.  The new simplex takes the next available
                    // index, and the canonical gluing becomes the identity.
                    image[dest] = nextUnusedSimp++;
                    preImage[image[dest]] = dest;
                    vertexMap[dest] = vertexMap[simpSrc] *
                        s->adjacentGluing(facetSrc).inverse();

                    facetAction[facetPos++] = 1;
                    continue;
                }

                // It's a simplex we've seen before.  Record the gluing.
                joinDest[joinPos] = image[dest];
                joinGluing[joinPos] = (vertexMap[dest] *
                    s->adjacentGluing(facetSrc) * vertexMap[simpSrc].inverse()).
                    index();
                ++joinPos;

                facetAction[facetPos++] = 2;
            }
        }

        // We have all we need.  Pack it all together into a string.
        // We need to encode:
        // - the number of simplices in this component;
        // - facetAction[i], 0 <= i < facetPos;
        // - joinDest[i], 0 <= i < joinPos;
        // - joinGluing[i], 0 <= i < joinPos.
        std::string ans;

        // Keep it simple for small triangulations (1 character per integer).
        // For large triangulations, start with a special marker followed by
        // the number of chars per integer.
        size_t nCompSimp = simpImg;
        unsigned nChars;
        if (nCompSimp < 63)
            nChars = 1;
        else {
            nChars = 0;
            size_t tmp = nCompSimp;
            while (tmp > 0) {
                tmp >>= 6;
                ++nChars;
            }

            ans = IsoSigHelper::SCHAR(63);
            ans += IsoSigHelper::SCHAR(nChars);
        }

        // Off we go.
        size_t i;
        IsoSigHelper::SAPPEND(ans, nCompSimp, nChars);
        for (i = 0; i < facetPos; i += 3)
            IsoSigHelper::SAPPENDTRITS(ans, facetAction + i,
                (facetPos >= i + 3 ? 3 : facetPos - i));
        for (i = 0; i < joinPos; ++i)
            IsoSigHelper::SAPPEND(ans, joinDest[i], nChars);
        for (i = 0; i < joinPos; ++i)
            IsoSigHelper::SAPPEND(ans, joinGluing[i],
                IsoSigHelper::CHARS_PER_PERM<dim>());

        // Record the canonical isomorphism if required.
        if (relabelling)
            for (i = 0; i < nCompSimp; ++i) {
                relabelling->simpImage(i) = image[i];
                relabelling->facetPerm(i) = vertexMap[i];
            }

        // Done!
        delete[] image;
        delete[] vertexMap;
        delete[] preImage;
        delete[] facetAction;
        delete[] joinDest;
        delete[] joinGluing;

        return ans;                
    }

public:
    //Need to debug segfault
    static std::string computeSignature(Triangulation<3>* triangulation) {
        std::vector<std::string> signatures;
        std::vector<SimplexInfo> properties;
        for (int i = 0; i < triangulation->size(); i++) {
            Simplex<3>* tetrahedra = triangulation->simplex(i);
            properties.emplace_back(SimplexInfo(tetrahedra, triangulation->size()));
        }
        std::sort(properties.begin(), properties.end());
        //Iterate through and partition into runs
        int prev = 0;
        int current = 1;
        
        std::vector<int> partition;
        for (int i = 1; i < properties.size(); i++) {
            if (properties[prev] == properties[i]) {
                current++;
            } else {
                partition.push_back(current);
                prev = i;
                current = 1;
            }
        }
        partition.push_back(current);
        //Use best partition (Partition requiring minimal tetrahedra)
        int index = 0;
        int bestIndex = 0; //Starting index for best tetrahedra
        int partitionIndex = 0; //Partition index for best tetrahedra
        int best = INT32_MAX;
        for (int i = 0; i < partition.size(); i++) {
            int count = 0;
            for (int j = 0; j < partition[i]; j++) {
                count += properties[index].numOrderings();
                index++;
            }
            if (count < best) {
                best = count;
                bestIndex = index - partition[i];
                partitionIndex = i;
            }
            best = std::min(count, best);
        }
        std::string ans;
        for (int i = 0; i < partition[partitionIndex]; i++) {
            auto perms = properties[bestIndex + i].getAllPerms();
            for (auto perm : perms) {
                std::string curr = isoSigFrom(triangulation, triangulation->simplex(bestIndex + i)->index(), 
                    Perm<4>::atIndex(perm), (Isomorphism<3>*) nullptr);
                if (ans.size() == 0) {
                    ans = curr;
                } else {
                    ans = std::min(ans, curr);
                }
            }    
        }
        return ans;
    }
};
#endif