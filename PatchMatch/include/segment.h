/*******************************************************************************
 * segments.h - segmentation tools
 *******************************************************************************
 * Add license here...
 *******************************/

#ifndef SEGMENTS_H
#define	SEGMENTS_H

#include <algebra.h>
#include <nnf.h>
#include <set>

namespace pm {

    template <typename T>
    inline bool highEnough(T t);

    template <> inline bool highEnough<int>(int i) {
        return i > 0;
    }

    template <> inline bool highEnough<float>(float f) {
        return f >= 0.5;
    }

    /**
     * \brief Minimal disjoint-set implementation for labeling
     *
     * \see http://en.wikipedia.org/wiki/Disjoint-set_data_structure
     */
    struct union_set {
        typedef union_set uset;
        // the data:
        uset *parent;
        int rank;

        union_set(){
            reset();
        }
        /// Actively get the root element
		
		void reset() {
			rank = 0;
			parent = this;
		}

        uset * root() {
            if (parent != this) parent = parent->root();
            return parent;
        }
        /// Merge two sets

        inline void link(union_set &s) {
            union_set *r1 = root();
            union_set *r2 = s.root();
            // if they are the same, nothing changes
            if (r1 == r2) return;

            // not in the same set => merge sets
            // /!\ do not use this->parent or s.parent as this / s may not be the roots!
            if (rank < s.rank) {
                r1->parent = r2;
            } else if (rank > s.rank) {
                r2->parent = r1;
            } else {
                r2->parent = r1;
                r1->rank += 1;
            }
        }
    };

    template <typename T>
    class Segmentation {
    public:
        typedef size_t Label;

        static const int ISOLATED = 0;
        static const int UP = 1 << 0;
        static const int DOWN = 1 << 1;
        static const int LEFT = 1 << 2;
        static const int RIGHT = 1 << 3;
        static const int INSIDE = UP | DOWN | LEFT | RIGHT;

        struct Node {
            int type;
            T data;

            Node() : type(ISOLATED), s() {
            }

            inline bool isBoundary() const {
                return (type & INSIDE) != INSIDE;
            }

            inline bool isInside() const {
                return (type & INSIDE) == INSIDE;
            }

            inline bool isUp() const {
                return (type & UP) == UP;
            }

            inline bool isLeft() const {
                return (type & LEFT) == LEFT;
            }

            inline bool isDown() const {
                return (type & DOWN) == DOWN;
            }

            inline bool isRight() const {
                return (type & RIGHT) == RIGHT;
            }

            inline Label label() {
                return reinterpret_cast<Label> (s.root());
            }

            inline void link(Node &p) {
                s.link(p.s);
            }
			
			inline void unlink() {
				s.reset();
			}

        private:
            union_set s;
            Node(const Node&);
        };

        const union {
			int width;
			int cols;
		};
		const union {
			int height;
			int rows;
		};

        Segmentation(int w, int h) : width(w), height(h), data(NULL) {
            data = new Node[h * w]();
        }
		
		void reset() {
			if(data == NULL) return;
			for(int i = 0, N = width * height; i < N; ++i) data[i].unlink();
		}

        ~Segmentation() {
            if (data != NULL) delete[] data;
        }

        inline const Node &at(int y, int x) const {
            return data[width * y + x];
        }

        inline Node &at(int y, int x) {
            return data[width * y + x];
        }

    private:
        Node *data;
    };
	
	struct NoData {};
	typedef Segmentation<NoData> BaseSegmentation;

    template <typename D, typename Patch, typename Scalar>
    Segmentation<D>* segment(NearestNeighborField<Patch, Scalar> *field, Segmentation<D>* segments = NULL) {
        typedef NearestNeighborField<Patch, Scalar> NNF;
        // segment and look at the shifts
        using patch::coherence;
        typedef Segmentation<D> Segments;
        typedef typename Segments::Node Node;

		// create object if needed
		if(segments == NULL){
			segments = new Segments(field->width, field->height);
		}
		
        // label the segments first
        for (int y = 0; y < field->height; ++y) {
            for (int x = 0; x < field->width; ++x) {
                Node &s = segments->at(y, x);
                const Patch &p = field->get(y, x);
                // by default, the pixel is isolated
                s.type = Segments::ISOLATED;
                if (y > 0) {
                    if (highEnough(coherence(p, field->get(y - 1, x), -1, 0))) {
                        s.type |= Segments::UP;
                        // update neighboring segment pixel
                        Node &n = segments->at(y - 1, x);
                        n.type |= Segments::DOWN;
                        // we link their label
                        s.link(n);
                    }
                }
                if (x > 0) {
                    if (highEnough(coherence(p, field->get(y, x - 1), 0, -1))) {
                        s.type |= Segments::LEFT;
                        // update neighboring segment pixel
                        Node &n = segments->at(y, x - 1); // neighboring segment pixel
                        n.type |= Segments::RIGHT;
                        // we link their label
                        s.link(n);
                    }
                }
            }
        }
        return segments;
    }
}

#endif	/* SEGMENTS_H */