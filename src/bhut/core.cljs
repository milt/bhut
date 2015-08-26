(ns ^:figwheel-always bhut.core
    (:require [reagent.core :as r]
              [reagent-forms.core :refer [bind-fields]]
              [clojure.zip :as zip]))

(enable-console-print!)


(defonce grav-constant (- 2000))
(defonce θ 0.5)
(defonce eta 3E4 ;; 10
  )

(defn square-points [min-x max-x
                     min-y max-y]
  (let [size-diff (- (.abs js/Math (- max-x min-x))
                     (.abs js/Math (- max-y min-y)))
        half-diff (/ size-diff 2)]
    (if (> size-diff 0)
      [min-x
       max-x
       (- min-y half-diff)
       (+ max-y half-diff)]
      [(+ min-x half-diff)
       (- max-x half-diff)
       min-y
       max-y])))


(defn add-bodies [x1 y1 m1
                  x2 y2 m2]
  (let [m (+ m1 m2)]
    [(/ (+ (* x1 m1)
           (* x2 m2)) m)
     (/ (+ (* y1 m1)
           (* y2 m2)) m)
     m]))

(defn calc-force [m1 m2 dx dy d]
  (let [f (/ (* grav-constant m1 m2) (+ (* d d) (* eta eta)))]
    [(/ (* f dx) d)
     (/ (* f dy) d)]))

{:x (- (rand-int 1000) 500)
 :y (- (rand-int 1000) 500)
 :vx 0.0
 :vy 0.0
 :m (+ 100 (rand-int 1000))}


(defrecord Body [x y vx vy m])


(declare new-bhtree tree-zipper smart-tree-zipper)

(defrecord BHTree [min-x
                   max-x
                   min-y
                   max-y
                   size
                   calc-size

                   com-x
                   com-y
                   mass

                   ^BHTree nw
                   ^BHTree ne
                   ^BHTree sw
                   ^BHTree se
                   bodies
                   forces
                   e-body]
  Object

  (children [this]
    (filter identity [nw ne sw se]))

  (internal? [this]
    (not e-body)
    #_(seq (.children this)))

  (contains-point? [_ {:keys [x y] :as body}]
    (and (<= min-x x)
         (>= max-x x)
         (<= min-y y)
         (>= max-y y)))

  (same-region? [this another-tree]
    (and (= min-x (:min-x another-tree))
         (= max-x (:max-x another-tree))
         (= min-y (:min-y another-tree))
         (= max-y (:max-y another-tree))))

  (sub-region [_ quad]
    (let [child-size (/ size 2)
          child-calc-size (* calc-size 2)
          [c-min-x c-max-x
           c-min-y c-max-y]
          (case quad
            :nw [min-x
                 (+ min-x child-size)
                 min-y
                 (+ min-y child-size)]
            :ne [(+ min-x child-size)
                 max-x
                 min-y
                 (+ min-y child-size)]
            :sw [min-x
                 (+ min-x child-size)
                 (+ min-y child-size)
                 max-y]
            :se [(+ min-x child-size)
                 max-x
                 (+ min-y child-size)
                 max-y])]
      (BHTree.
       c-min-x
       c-max-x
       c-min-y
       c-max-y

       child-size
       child-calc-size

       nil nil nil
       nil nil nil nil
       nil
       nil
       nil)))

  (branch [this {:keys [x y] :as body}]
    (some
     (fn [k]
       (let [quad (or (k this) (.sub-region this k))]
         (when (.contains-point? quad body)
           [k quad])))
     [:nw :ne :sw :se]))

  (insert [this {:keys [x y m] :as body}]
    (cond
      ;; if it's not in this region, we don't care!
      (not (.contains-point? this body))
      this

      ;; if this is a fresh node, it has no mass, update it!
      (nil? mass)
      (assoc this
             :com-x x
             :com-y y
             :mass m
             :e-body body)

      ;; if this is an internal node, add c-o-m find the right sub and insert
      (.internal? this)
      (let [[k tree] (.branch this body)
            [new-com-x
             new-com-y
             new-mass] (add-bodies x y m
                                   com-x com-y mass)]
        (assoc this
               :com-x new-com-x
               :com-y new-com-y
               :mass new-mass
               k (.insert tree body)))

      ;; Otherwise, this is an external node and we need to add c-o-m, then propagate both
      :else
      (let [[ak atree] (.branch this {:x com-x :y com-y})
            [bk btree] (.branch this body)

            [new-com-x
             new-com-y
             new-mass] (add-bodies x y m
                                   com-x com-y mass)]
        (cond-> (assoc this
                       :com-x new-com-x
                       :com-y new-com-y
                       :mass new-mass
                       :e-body nil)
          (= ak bk) (assoc ak (-> atree
                                  (.insert {:x com-x :y com-y :m mass})
                                  (.insert body)))
          (not= ak bk) (assoc
                        ak (.insert atree {:x com-x :y com-y :m mass})
                        bk (.insert btree body))))))

  (bh-tree-seq [this]
    (tree-seq #(.internal? %)
              #(.children %)
              this))

  (load-bodies [this]
    (reduce
     (fn [t b]
       (.insert t b))
     this
     bodies))

  (clear-tree [this]
    (assoc this
           :com-x nil
           :com-y nil
           :mass nil
           :e-body nil
           :nw nil
           :ne nil
           :sw nil
           :se nil))

  (forces-seq [this {:keys [x y m] :as body} theta]
    (let [calc-d (fn [{:keys [com-x com-y mass size] :as node}]
                   (let [dx (- x com-x)
                         dy (- y com-y)
                         d (.sqrt js/Math (+ (* dx dx) (* dy dy)))]
                     [dx dy d]))
          ;; branch? (fn [{:keys [com-x com-y mass size e-body] :as node}]
          ;;           (and (not e-body)
          ;;                (let [[dx dy d] (calc-d node)]
          ;;                  (> (/ size d) theta))))
          ;; children (fn [node]
          ;;            (.children node))
          walk (fn walk [node]
                 (lazy-seq
                  (let [{:keys [mass size e-body]} node
                        [dx dy d] (calc-d node)
                        branch? (and (not e-body)
                                     (> (/ size d) theta))]

                    (cond
                      branch?
                      (mapcat walk (.children node))

                      (= d 0)
                      nil

                      :else
                      (cons
                       (calc-force mass m dx dy d)
                       (when branch?
                         (mapcat walk (.children node))))
                      ))))]
      (if (.contains-point? this body)
        (doall (walk this))
        (list
         (let [[dx dy d] (calc-d this)]
           (calc-force m mass dx dy d))))))

  (get-node-forces [this theta]
    (into {}
          (for [{:keys [x y m] :as body} bodies]
            [body
             (reduce
              (fn [[fx fy] [dfx dfy]]
                [(+ fx dfx)
                 (+ fy dfy)])
              [0 0]
              (.forces-seq this body theta))])))

  (get-all-forces [this theta]
    (-> this
        .clear-tree
        .load-bodies
        (assoc :forces (.get-node-forces this theta))))

  (update-velocities [this dt]
    (if (= forces {})
      this
      (assoc this
            :bodies
            (into []
                  (for [{:keys [x y vx vy m] :as body} bodies
                        :let [[fx fy] (get forces body)]]
                    (assoc
                     body
                     :vx (+ vx (/ (* dt fx) m))
                     :vy (+ vy (/ (* dt fy) m)))))
            :forces {})))

  (update-positions [this dt]
    (assoc this
           :bodies
           (into []
                 (for [{:keys [x y vx vy m] :as body} bodies]
                   (assoc
                    body
                    :x (+ x (* vx dt))
                    :y (+ y (* vy dt)))))))

  (add-body [this {:keys [x y vx vy m] :as body}]
    (assoc
     (.insert this body)
     :bodies (conj bodies body)))

  (add-bodies [this bodies]
    (reduce
     (fn [t b]
       (.add-body t b))
     this
     bodies)))


(defn new-bhtree [min-x max-x
                  min-y max-y
                  bodies]
  (let [
        ;; [min-x max-x
        ;;  min-y max-y] (square-points min-x max-x
        ;;                              min-y max-y)
         size (.abs js/Math (- max-x min-x))
         calc-size (/ 1 size)]
    (BHTree.
     min-x
     max-x
     min-y
     max-y
     size
     calc-size
     nil nil nil ;; body info
     nil nil nil nil ;; quadrants
     bodies
     {} ;; forces
     nil ;; external-body
     )))


(defn rand-points [n]
  (for [i (range n)]
    (map->Body
     {:x (- (rand-int 1000) 500)
      :y (- (rand-int 1000) 500)
      :vx 0.0
      :vy 0.0
      :m (+ 100 (rand-int 1000))})))

(defn update-pos [tree dt theta]
  (-> tree
      (.get-all-forces theta)
      (.update-velocities dt)
      (.update-positions dt)
      ))


(defn upd-f [tree dt theta]
  (.get-all-forces tree theta))
(defn upd-v [tree dt theta]
  (.update-velocities tree dt))
(defn upd-p [tree dt theta]
  (.update-positions tree dt))


(defn add-b [tree b]
  (.add-body tree b))

(defn add-bs [tree bs]
  (.add-bodies tree bs))

(def knob-template
  [:div
   [:label {:for "dt"} "dt"]
   [:input {:field :numeric :id :dt}]
   [:label {:for "theta"} "theta"]
   [:input {:field :numeric :id :theta}]
   [:label {:for "debug"} "debug"]
   [:input {:field :checkbox :id :debug}]])

(defn- cycle-fns [fn-atom]
  (concat (rest fn-atom)
          (list (first fn-atom))))

(defn bhut []
  (let [knobs (r/atom {:dt 0.5
                     :theta 0.5
                     :debug false})
        tree (r/atom (new-bhtree (- 500) 500 (- 500) 500 []))
        points (r/cursor tree [:bodies])
        running? (r/atom false)
        step-fn #(let [{:keys [dt theta]} @knobs]
                   (swap! tree update-pos dt theta))]
    (js/setInterval
     (fn []
       (when @running?
         (r/next-tick step-fn))) 60)
    (fn []
      (let [
            ;; tree (.load-bodies (new-bhtree (- 500) 500 (- 500) 500) @points)
            ;; step-fn
            ]
        ;; (print (count @acting-nodes))
        [:div
         [bind-fields knob-template knobs]
         [:button {:on-click (fn [e]
                               (swap! tree add-b (first (rand-points 1))))}
          "one"]
         [:button {:on-click (fn [e]
                               (swap! tree add-bs  (rand-points 100)))}
          "100"]
         [:button {:on-click #(time (step-fn))}
          "step"]
         [:button {:on-click #(reset! tree (new-bhtree (- 500) 500 (- 500) 500 []))}
          "clear"]
         [:button {:on-click #(swap! running? not)}
          (if @running?
            "Stop"
            "Start")]
         [:button {:on-click #(time (-> @tree
                                         .clear-tree
                                         .load-bodies
                                         (.get-node-forces (:theta @knobs))))}
          "dump node forces"]
         [:br]
         [:svg
          {:width 1000 :height 1000
           :on-double-click (fn [e]
                       (let [
                             t (.-currentTarget e)
                             dim (.getBoundingClientRect t)
                             svg-x (- (.-clientX e) (.-left dim))
                             svg-y (- (.-clientY e) (.-top dim))]
                         (swap! tree add-b (map->Body
                                              {:x (- svg-x 500)
                                               :y (- svg-y 500)
                                               :vx 0.0
                                               :vy 0.0
                                               :m 1000}))))}
          [:g
           {:transform "translate(500,500)"}
           (when (:debug @knobs)
             (into [:g.tree]
                  (for [{:keys [min-x min-y
                                size
                                com-x com-y mass]} (.bh-tree-seq @tree)]
                    [:g.node
                     [:circle.center-of-mass
                      {:cx com-x
                       :cy com-y
                       :r (/ mass 100)
                       :style
                       {:stroke "blue"
                        :stroke-width 2
                        :fill-opacity 0.0}}]
                     [:rect.quad
                      {:x min-x
                       :y min-y
                       :width size
                       :height size
                       :style {:stroke "red"
                               :stroke-width 5
                               :fill-opacity 0.0}}]])))
           (into [:g.points]
                 (for [{:keys [x y m] :as body} @points]
                   [:circle
                    {
                     :on-click
                     #(-> @tree
                         .clear-tree
                         .load-bodies
                         (.forces-seq body (:theta @knobs)))
                     :cx x
                     :cy y
                     :r (/ m 100)
                     :fill "green"}]))]]]))))








;; define your app data so that it doesn't get over-written on reload

(defonce app-state (atom {:text "Hello world!"}))


(defn mount []
  (r/render-component [bhut]
                      (.-body js/document)))

(mount)

(defn on-js-reload []
  ;; optionally touch your app-state to force rerendering depending on
  ;; your application
  ;; (swap! app-state update-in [:__figwheel_counter] inc)
  (mount)
)
