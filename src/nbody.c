#include "nbody.h"
#include "common.h"
#include "sources.h"
#include "utils.h"

#include <math.h>

// NBODY_G: gravitational constant, hand-tuned with the masses in
// init_sources() for noticeable but non-chaotic motion.
// NBODY_SOFTENING: minimum effective distance between two sources, keeps
// force (proportional to 1/dist^2) from blowing up on a close encounter.
// NBODY_VEL_MAX: max velocity per axis (cells/frame), caps the slingshot
// effect from a close encounter so the integration stays stable.
// NBODY_RESTITUTION: fraction of velocity kept on a wall bounce (below 1 =
// inelastic).
#define NBODY_G            0.3f
#define NBODY_SOFTENING    6.0f
#define NBODY_VEL_MAX      5.0f
#define NBODY_RESTITUTION  0.9f

// Barnes-Hut: max angular size (node width / distance) before a node is
// approximated as a single point mass instead of descending into its
// children. Smaller = more accurate (closer to true pairwise gravity) but
// slower.
#define BH_THETA      0.5f
#define BH_MAX_NODES  (SOURCES_MAX * 8)
#define BH_MAX_DEPTH  30

// One quadtree node: a square region, its total mass and center of mass,
// and up to 4 children (one per quadrant). A leaf either holds exactly one
// source (`body` >= 0) or is empty (`body` == -1, mass == 0).
typedef struct {
    float cx, cy, half;
    float mass, com_x, com_y;
    int   body;
    int   children[4];
} QuadNode;

static QuadNode bh_pool[BH_MAX_NODES];
static int      bh_count;

static int bh_new_node(float cx, float cy, float half)
{
    QuadNode *n = &bh_pool[bh_count];
    n->cx = cx; n->cy = cy; n->half = half;
    n->mass = 0.0f; n->com_x = 0.0f; n->com_y = 0.0f;
    n->body = -1;
    n->children[0] = n->children[1] = n->children[2] = n->children[3] = -1;
    return bh_count++;
}

// Quadrant index (0-3) of point (x, y) relative to this node's center.
static int bh_quadrant(const QuadNode *n, float x, float y)
{
    int q = 0;
    if (x >= n->cx) q |= 1;
    if (y >= n->cy) q |= 2;
    return q;
}

static void bh_subdivide(int idx)
{
    QuadNode *n = &bh_pool[idx];
    float h = n->half * 0.5f;
    n->children[0] = bh_new_node(n->cx - h, n->cy - h, h);
    n->children[1] = bh_new_node(n->cx + h, n->cy - h, h);
    n->children[2] = bh_new_node(n->cx - h, n->cy + h, h);
    n->children[3] = bh_new_node(n->cx + h, n->cy + h, h);
}

// Inserts one source into the tree, splitting leaves into 4 children as
// needed, updating mass/center of mass on every visited node on the way
// back up.
static void bh_insert(int idx, int body, InkSource *sources, int depth)
{
    QuadNode *n = &bh_pool[idx];
    float bx = sources[body].pos_x, by = sources[body].pos_y;
    float bm = sources[body].mass;

    if (n->body == -1 && n->children[0] == -1 && n->mass == 0.0f) {
        n->body = body;
        n->mass = bm;
        n->com_x = bx;
        n->com_y = by;
        return;
    }

    if (n->children[0] == -1) {
        // Leaf already holding one source: split, then push both down.
        // Depth guard: nearly-coincident sources would otherwise recurse
        // forever, so past BH_MAX_DEPTH just merge into this leaf's
        // combined mass/center of mass.
        if (depth >= BH_MAX_DEPTH) {
            float total = n->mass + bm;
            n->com_x = (n->com_x * n->mass + bx * bm) / total;
            n->com_y = (n->com_y * n->mass + by * bm) / total;
            n->mass = total;
            return;
        }
        int old_body = n->body;
        n->body = -1;
        bh_subdivide(idx);
        bh_insert(n->children[bh_quadrant(n, sources[old_body].pos_x,
                                          sources[old_body].pos_y)],
                 old_body, sources, depth + 1);
        bh_insert(n->children[bh_quadrant(n, bx, by)], body, sources, depth + 1);
    } else {
        bh_insert(n->children[bh_quadrant(n, bx, by)], body, sources, depth + 1);
    }

    // Accumulate this source's mass/position into the node's combined
    // mass/center of mass.
    float total = n->mass + bm;
    n->com_x = (n->com_x * n->mass + bx * bm) / total;
    n->com_y = (n->com_y * n->mass + by * bm) / total;
    n->mass = total;
}

// Walks the tree accumulating gravitational acceleration on `body`. A node
// is treated as a single point mass once it's a leaf, or once its angular
// size is below BH_THETA; otherwise recurses into its 4 children.
static void bh_accumulate(int idx, int body, const InkSource *sources,
                          float *acc_x, float *acc_y)
{
    const QuadNode *n = &bh_pool[idx];
    // Empty node, or the node's only body is the one we're computing for
    // (no self-attraction).
    if (n->mass <= 0.0f || n->body == body) return;

    float dx = n->com_x - sources[body].pos_x;
    float dy = n->com_y - sources[body].pos_y;
    // Squared distance to the node's center of mass, softened by
    // NBODY_SOFTENING^2 so it never approaches 0 (which would send the
    // force below toward infinity as two bodies converge).
    float dist2 = dx * dx + dy * dy + NBODY_SOFTENING * NBODY_SOFTENING;

    // Barnes-Hut opening criterion: node width / distance < BH_THETA. A
    // leaf always qualifies (no children to open). Passing this treats the
    // node as one point mass; failing it opens the node into its children.
    int is_leaf = (n->children[0] == -1);
    if (is_leaf || (n->half * 2.0f) / sqrtf(dist2) < BH_THETA) {
        // Newtonian gravity, F = G*m1*m2/d^2, directed from body toward
        // the node's center of mass. n->mass already sums every body
        // inside the node, so this one calculation stands in for all of
        // them. Multiplying by dx/dy (not a unit vector) folds direction
        // and magnitude into the same computation.
        float inv_dist = 1.0f / sqrtf(dist2);
        float factor = NBODY_G * n->mass * inv_dist * inv_dist * inv_dist;
        *acc_x += factor * dx;
        *acc_y += factor * dy;
        return;
    }

    // Node failed the opening criterion: recurse into its 4 quadrants.
    for (int c = 0; c < 4; c++) {
        if (n->children[c] != -1) {
            bh_accumulate(n->children[c], body, sources, acc_x, acc_y);
        }
    }
}

// Advances one frame of the n-body system. Each source attracts every
// other source (F = G*m1*m2/d^2, softened) and bounces off the grid's
// edges to stay on screen.
//
// Gravity is computed via a Barnes-Hut quadtree instead of direct pairwise
// summation: O(count log count) instead of O(count^2). The tree is rebuilt
// from scratch every frame since sources move every frame.
//
// Runs on a per-frame time unit, not the fluid solver's dt: source movement
// is a separate system that only happens to share the grid.
void update_nbody_sources(InkSource *sources, int count, int resolution)
{
    float accel_x[SOURCES_MAX] = { 0 };
    float accel_y[SOURCES_MAX] = { 0 };
    int   i;
    const int radius = source_radius(resolution);
    // Sources can't leave the injection neighborhood (source_radius() in
    // sources.h) without losing valid interior cells.
    const float lower_bound = (float)(radius + 1);
    const float upper_bound = (float)(resolution - radius);

    // Build the tree. Root covers the whole grid.
    bh_count = 0;
    int root = bh_new_node((float)resolution * 0.5f, (float)resolution * 0.5f,
                           (float)resolution * 0.5f);
    for (i = 0; i < count; i++) {
        bh_insert(root, i, sources, 0);
    }

    for (i = 0; i < count; i++) {
        bh_accumulate(root, i, sources, &accel_x[i], &accel_y[i]);
    }

    for (i = 0; i < count; i++) {
        sources[i].vel_x = clamp(sources[i].vel_x + accel_x[i],
                                 -NBODY_VEL_MAX, NBODY_VEL_MAX);
        sources[i].vel_y = clamp(sources[i].vel_y + accel_y[i],
                                 -NBODY_VEL_MAX, NBODY_VEL_MAX);

        sources[i].pos_x += sources[i].vel_x;
        sources[i].pos_y += sources[i].vel_y;

        // Damped bounce off the grid's edges, keeps the source on screen.
        if (sources[i].pos_x < lower_bound) {
            sources[i].pos_x = lower_bound;
            sources[i].vel_x = -sources[i].vel_x * NBODY_RESTITUTION;
        } else if (sources[i].pos_x > upper_bound) {
            sources[i].pos_x = upper_bound;
            sources[i].vel_x = -sources[i].vel_x * NBODY_RESTITUTION;
        }

        if (sources[i].pos_y < lower_bound) {
            sources[i].pos_y = lower_bound;
            sources[i].vel_y = -sources[i].vel_y * NBODY_RESTITUTION;
        } else if (sources[i].pos_y > upper_bound) {
            sources[i].pos_y = upper_bound;
            sources[i].vel_y = -sources[i].vel_y * NBODY_RESTITUTION;
        }
    }
}
